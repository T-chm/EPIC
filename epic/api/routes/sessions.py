"""Session management endpoints."""

from __future__ import annotations

from dataclasses import asdict

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from epic.api.dependencies import (
    create_session,
    delete_session,
    get_session,
    list_sessions,
    switch_session_provider,
)

router = APIRouter()


class CreateSessionRequest(BaseModel):
    provider: str = "ollama"
    model: str | None = None


class SwitchProviderRequest(BaseModel):
    provider: str
    model: str | None = None


@router.post("", status_code=201)
async def create_new_session(req: CreateSessionRequest) -> dict:
    """Create a new chat session."""
    try:
        info = create_session(req.provider, req.model)
        return asdict(info)
    except (ValueError, ImportError) as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("")
async def list_all_sessions() -> list[dict]:
    """List all active sessions."""
    return [asdict(s) for s in list_sessions()]


@router.get("/{session_id}")
async def get_session_info(session_id: str) -> dict:
    """Get session details with message history."""
    session = get_session(session_id)
    if not session:
        raise HTTPException(status_code=404, detail="Session not found")

    info = asdict(session["info"])
    engine = session["engine"]
    info["messages"] = [
        {"role": m.role.value, "content": m.content} for m in engine.history
    ]
    return info


@router.delete("/{session_id}", status_code=204)
async def delete_existing_session(session_id: str) -> None:
    """Delete a session."""
    if not delete_session(session_id):
        raise HTTPException(status_code=404, detail="Session not found")


@router.patch("/{session_id}/provider")
async def switch_provider(session_id: str, req: SwitchProviderRequest) -> dict:
    """Switch the LLM provider for an existing session."""
    try:
        info = switch_session_provider(session_id, req.provider, req.model)
        if not info:
            raise HTTPException(status_code=404, detail="Session not found")
        return asdict(info)
    except (ValueError, ImportError) as e:
        raise HTTPException(status_code=400, detail=str(e))
