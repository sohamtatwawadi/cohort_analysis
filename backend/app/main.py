"""FastAPI application entry point."""
from __future__ import annotations

import logging
import re

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, Response
from fastapi.staticfiles import StaticFiles

from . import config, db
from .api.research_routes import router as research_router
from .api.routes import router
from .ingest import store, synthetic
from .research import jobs as research_jobs
from .services import library

log = logging.getLogger("germline")

app = FastAPI(
    title="ImpactOmics · Germline Cohort Analytics",
    description=("Subject-level germline cohort analytics. Every number is computed "
                 "server-side over a rebuildable analytics store derived from VariMAT."),
    version="1.0.0",
)

app.add_middleware(
    CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

app.include_router(router, prefix="/api")
app.include_router(research_router, prefix="/api/research")


@app.on_event("startup")
def startup() -> None:
    store.init()
    if store.is_empty():
        log.warning("Analytics store is empty — seeding synthetic data. "
                    "Run `python -m backend.cli ingest <dir>` to load VariMAT files.")
        synthetic.rebuild()
    library.seed_builtins()

    # Job state is in the database; the thread advancing it is not. A job that
    # was running when the process stopped has nobody left to finish it, and is
    # polled forever from the UI. Startup is the one moment nothing can be
    # genuinely in flight, so it is where that gets reconciled.
    interrupted = research_jobs.reconcile_interrupted()
    if interrupted:
        log.warning("Marked %d interrupted analysis job(s) as failed.", interrupted)


@app.get("/health")
def health():
    return {"ok": True, "store": db.table_counts()}


if config.FRONTEND_DIR.exists():
    @app.middleware("http")
    async def no_store_static(request, call_next):
        """Browsers cache ES modules aggressively, so an edited view keeps
        rendering the previous build until a manual hard reload. Costs nothing
        here — the assets are served from the same box."""
        response = await call_next(request)
        if request.url.path.startswith("/static") or request.url.path == "/":
            response.headers["Cache-Control"] = "no-store, must-revalidate"
        return response

    def _asset_version() -> str:
        """Newest mtime across the frontend sources.

        ES modules are cached by URL, and a `no-store` header does not evict an
        entry the browser already holds — which means an edited view can keep
        rendering the previous build indefinitely. Stamping the URL makes the
        cache key change whenever any asset changes, so there is nothing to
        evict.
        """
        newest = 0.0
        for p in config.FRONTEND_DIR.rglob("*"):
            if p.suffix in (".js", ".css"):
                newest = max(newest, p.stat().st_mtime)
        return str(int(newest))

    _IMPORT_RE = re.compile(r"(from\s+['\"])(\.{1,2}/[^'\"]+?\.js)(['\"])")

    @app.get("/static/js/{name}.js")
    def module(name: str):
        """Serve a frontend module with its sibling imports version-stamped.

        Stamping only the entry point is not enough: `app.js?v=N` still resolves
        `./views.js` to an unversioned URL, so an edited view keeps loading from
        cache while the shell reloads. Rewriting the specifiers propagates the
        version through the whole module graph.
        """
        path = config.FRONTEND_DIR / "js" / (name + ".js")
        if not path.is_file():
            raise HTTPException(404, "no such module")
        v = _asset_version()
        src = _IMPORT_RE.sub(lambda m: m.group(1) + m.group(2) + "?v=" + v + m.group(3),
                             path.read_text())
        return Response(src, media_type="application/javascript")

    app.mount("/static", StaticFiles(directory=str(config.FRONTEND_DIR)), name="static")

    @app.get("/")
    def index():
        html = (config.FRONTEND_DIR / "index.html").read_text()
        v = _asset_version()
        # Stamp every local asset, not one named entry point — the entry
        # module's name is a detail of the frontend, and hard-coding it means a
        # rename silently disables cache-busting.
        html = re.sub(r'(/static/[A-Za-z0-9_./-]+\.(?:css|js))',
                      lambda m: m.group(1) + "?v=" + v, html)
        return HTMLResponse(html)
