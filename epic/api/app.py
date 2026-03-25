"""FastAPI application factory."""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from epic.api.middleware import setup_middleware
from epic.config import get_settings

logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan: startup and shutdown logic."""
    settings = get_settings()
    logging.basicConfig(level=getattr(logging, settings.log_level.upper(), logging.INFO))
    logger.info("EPIC API starting up")

    # Pre-warm Ollama model to eliminate cold-start latency
    if settings.default_provider == "ollama":
        try:
            from epic.providers.ollama_provider import _get_client

            client = _get_client(settings.ollama_host)
            client.chat(
                model=settings.ollama_default_model,
                messages=[{"role": "user", "content": "hi"}],
                keep_alive=settings.ollama_keep_alive,
            )
            logger.info("Pre-warmed Ollama model: %s", settings.ollama_default_model)
        except Exception as e:
            logger.warning("Ollama pre-warm failed: %s", e)

    yield
    logger.info("EPIC API shutting down")


def create_app() -> FastAPI:
    """Create and configure the FastAPI application."""
    app = FastAPI(
        title="EPIC API",
        description="Engineering Probes via Instructing a Chatbot — REST + WebSocket API",
        version="2.0.0",
        lifespan=lifespan,
    )

    setup_middleware(app)

    # Register route modules
    from epic.api.routes import chat, code, health, providers, sessions

    app.include_router(health.router, prefix="/api", tags=["health"])
    app.include_router(providers.router, prefix="/api/providers", tags=["providers"])
    app.include_router(sessions.router, prefix="/api/sessions", tags=["sessions"])
    app.include_router(code.router, prefix="/api/code", tags=["code"])
    app.include_router(chat.router, tags=["chat"])

    # Serve React frontend build in production
    frontend_build = Path(__file__).parent.parent.parent / "frontend" / "dist"
    if frontend_build.exists():
        app.mount(
            "/",
            StaticFiles(directory=str(frontend_build), html=True),
            name="frontend",
        )
        logger.info("Serving frontend from %s", frontend_build)

    return app
