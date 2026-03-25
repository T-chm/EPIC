"""Provider management endpoints."""

from __future__ import annotations

from fastapi import APIRouter, HTTPException

from epic.config import get_settings
from epic.providers import get_provider

router = APIRouter()


@router.get("")
async def list_providers() -> list[dict]:
    """List available providers and their models."""
    settings = get_settings()
    results = []

    for name in ("openai", "anthropic", "ollama"):
        try:
            provider = get_provider(name, settings)
            models = provider.list_models()
            results.append({
                "name": name,
                "available": True,
                "models": models,
                "default_model": provider.model_name,
            })
        except (ValueError, ImportError):
            results.append({
                "name": name,
                "available": False,
                "models": [],
                "default_model": None,
            })

    return results


@router.get("/{provider_name}/models")
async def list_models(provider_name: str) -> dict:
    """List models for a specific provider."""
    settings = get_settings()
    try:
        provider = get_provider(provider_name, settings)
        models = provider.list_models()
        return {"provider": provider_name, "models": models}
    except (ValueError, ImportError) as e:
        raise HTTPException(status_code=404, detail=str(e))
