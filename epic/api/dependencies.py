"""Shared FastAPI dependencies: session and engine management."""

from __future__ import annotations

import uuid
from datetime import datetime, timezone

from epic.chat import ChatEngine
from epic.config import get_settings
from epic.models import SessionInfo
from epic.providers import get_provider

# In-memory session store (for production, swap to Redis/DB)
_sessions: dict[str, dict] = {}


def create_session(provider_name: str, model: str | None = None) -> SessionInfo:
    """Create a new chat session."""
    settings = get_settings()
    provider = get_provider(provider_name, settings, model=model)
    session_id = uuid.uuid4().hex[:12]

    engine = ChatEngine(
        provider=provider,
        auto_execute=settings.auto_execute_code,
        temperature=settings.default_temperature,
        max_tokens=settings.default_max_tokens,
        compression_threshold=settings.memory_compression_threshold,
        compression_keep_recent=settings.memory_keep_recent,
    )

    info = SessionInfo(
        session_id=session_id,
        provider_name=provider.name,
        model=provider.model_name,
        created_at=datetime.now(timezone.utc).isoformat(),
    )

    _sessions[session_id] = {
        "info": info,
        "engine": engine,
    }

    return info


def get_session(session_id: str) -> dict | None:
    """Get a session by ID. Returns None if not found."""
    return _sessions.get(session_id)


def get_engine(session_id: str) -> ChatEngine | None:
    """Get the ChatEngine for a session."""
    session = _sessions.get(session_id)
    return session["engine"] if session else None


def delete_session(session_id: str) -> bool:
    """Delete a session. Returns True if it existed."""
    return _sessions.pop(session_id, None) is not None


def list_sessions() -> list[SessionInfo]:
    """List all active sessions."""
    result = []
    for s in _sessions.values():
        info: SessionInfo = s["info"]
        info.message_count = len(s["engine"].history)
        result.append(info)
    return result


def switch_session_provider(
    session_id: str, provider_name: str, model: str | None = None
) -> SessionInfo | None:
    """Switch the provider for an existing session."""
    session = _sessions.get(session_id)
    if not session:
        return None

    settings = get_settings()
    new_provider = get_provider(provider_name, settings, model=model)
    engine: ChatEngine = session["engine"]
    engine.switch_provider(new_provider)

    info: SessionInfo = session["info"]
    info.provider_name = new_provider.name
    info.model = new_provider.model_name
    return info
