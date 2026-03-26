"""FastAPI application for ALASCCA-classify web interface."""

from __future__ import annotations

import uuid
from collections import OrderedDict
from pathlib import Path

from fastapi import APIRouter, FastAPI, File, Form, Request, UploadFile
from fastapi.responses import HTMLResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from .. import __version__
from ..classify import classify_sample
from ..io import parse_maf, parse_seg
from ..models import SampleInput, SampleMetadata
from ..report import result_to_dict, result_to_pdf

WEB_DIR = Path(__file__).parent
TEMPLATE_DIR = WEB_DIR / "templates"
STATIC_DIR = WEB_DIR / "static"

MAX_RESULTS = 500


def _find_demo_dir() -> Path:
    candidates = [
        Path(__file__).resolve().parent.parent.parent.parent / "demo",
        Path("/app/demo"),
        Path.cwd() / "demo",
    ]
    for p in candidates:
        if p.is_dir():
            return p
    return candidates[0]


DEMO_DIR = _find_demo_dir()


class ResultStore:
    """Simple in-memory result store with FIFO eviction."""

    def __init__(self, maxsize: int = MAX_RESULTS):
        self._store: OrderedDict[str, dict] = OrderedDict()
        self._maxsize = maxsize

    def put(self, data: dict) -> str:
        rid = uuid.uuid4().hex[:12]
        if len(self._store) >= self._maxsize:
            self._store.popitem(last=False)
        self._store[rid] = data
        return rid

    def get(self, rid: str) -> dict | None:
        return self._store.get(rid)


# --- Routes ---

classify_router = APIRouter()
methods_router = APIRouter()


@classify_router.get("/", response_class=HTMLResponse)
async def index(request: Request):
    templates = request.app.state.templates
    return templates.TemplateResponse(request, "classify.html")


@classify_router.get("/health")
async def health():
    return {"status": "ok", "version": __version__}


@classify_router.post("/classify", response_class=HTMLResponse)
async def classify_route(
    request: Request,
    maf_file: UploadFile = File(...),
    seg_file: UploadFile | None = File(None),
    sample_id: str = Form(""),
    ptnm_stage: str = Form(""),
    tumor_location: str = Form(""),
    msi_status: str = Form(""),
    tumor_purity: str = Form(""),
):
    templates = request.app.state.templates
    results_store = request.app.state.results

    try:
        maf_content = (await maf_file.read()).decode("utf-8")
        variants = parse_maf(maf_content)

        # Build metadata
        sid = sample_id.strip() or _infer_sample_id(maf_content, maf_file.filename)
        metadata = SampleMetadata(
            sample_id=sid,
            ptnm_stage=ptnm_stage.strip(),
            tumor_location=tumor_location.strip(),
            msi_status=msi_status.strip(),
            tumor_purity=_parse_float(tumor_purity),
        )

        # Parse SEG if provided
        segments = []
        if seg_file and seg_file.filename:
            seg_content = (await seg_file.read()).decode("utf-8")
            segments = parse_seg(seg_content)

        sample = SampleInput(metadata=metadata, variants=variants, segments=segments)
        result = classify_sample(sample)

        result_data = result_to_dict(result)
        rid = results_store.put(result_data)

        return templates.TemplateResponse(request, "result.html", {
            "r": result_data,
            "result_id": rid,
        })
    except Exception as e:
        return templates.TemplateResponse(request, "error.html", {
            "error": str(e),
        })


@classify_router.get("/demo/{sample_name}", response_class=HTMLResponse)
async def demo(request: Request, sample_name: str):
    templates = request.app.state.templates
    results_store = request.app.state.results

    maf_path = DEMO_DIR / f"{sample_name}.maf"
    seg_path = DEMO_DIR / f"{sample_name}.seg"

    if not maf_path.exists():
        return templates.TemplateResponse(request, "error.html", {
            "error": f"Demo sample '{sample_name}' not found.",
        })

    try:
        variants = parse_maf(maf_path)

        segments = []
        if seg_path.exists():
            segments = parse_seg(seg_path)

        metadata = SampleMetadata(sample_id=sample_name)
        sample = SampleInput(metadata=metadata, variants=variants, segments=segments)
        result = classify_sample(sample)

        result_data = result_to_dict(result)
        rid = results_store.put(result_data)

        return templates.TemplateResponse(request, "result.html", {
            "r": result_data,
            "result_id": rid,
            "is_demo": True,
        })
    except Exception as e:
        return templates.TemplateResponse(request, "error.html", {
            "error": str(e),
        })


@classify_router.get("/results/{result_id}/pdf")
async def download_pdf(request: Request, result_id: str):
    result_data = request.app.state.results.get(result_id)
    if result_data is None:
        return HTMLResponse("Result not found or expired.", status_code=404)

    from ..models import ClassificationResult

    result = ClassificationResult(**result_data)
    pdf_bytes = result_to_pdf(result)
    sid = result_data.get("sample_id", "sample")
    return StreamingResponse(
        iter([pdf_bytes]),
        media_type="application/pdf",
        headers={"Content-Disposition": f'attachment; filename="{sid}_alascca_report.pdf"'},
    )


@classify_router.get("/results/{result_id}/json")
async def result_json(request: Request, result_id: str):
    result_data = request.app.state.results.get(result_id)
    if result_data is None:
        return HTMLResponse("Result not found or expired.", status_code=404)
    return result_data


@methods_router.get("/methods", response_class=HTMLResponse)
async def methods(request: Request):
    templates = request.app.state.templates
    return templates.TemplateResponse(request, "methods.html")


@methods_router.get("/install", response_class=HTMLResponse)
async def install(request: Request):
    templates = request.app.state.templates
    return templates.TemplateResponse(request, "install.html")


# --- App factory ---

def create_app() -> FastAPI:
    app = FastAPI(
        title="alascca-classify",
        description="PI3K pathway alteration classification for colorectal cancer",
        version=__version__,
    )

    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

    templates = Jinja2Templates(directory=str(TEMPLATE_DIR))
    templates.env.globals["version"] = __version__

    app.state.templates = templates
    app.state.results = ResultStore()

    app.include_router(classify_router)
    app.include_router(methods_router)

    return app


app = create_app()


def _parse_float(val: str) -> float | None:
    if not val or not val.strip():
        return None
    try:
        return float(val.strip())
    except ValueError:
        return None


def _infer_sample_id(maf_content: str, filename: str | None) -> str:
    if filename:
        return Path(filename).stem
    return "unknown"
