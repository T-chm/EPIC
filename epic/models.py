"""Shared data models for EPIC."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class Role(str, Enum):
    SYSTEM = "system"
    USER = "user"
    ASSISTANT = "assistant"


@dataclass
class Message:
    role: Role
    content: str


@dataclass
class StreamEvent:
    """Events yielded during streaming from ChatEngine."""

    event_type: str  # "token", "code_detected", "code_executing", "code_result", "done", "error"
    data: str

    def to_dict(self) -> dict[str, Any]:
        return {"type": self.event_type, "data": self.data}


@dataclass
class CodeResult:
    code: str
    stdout: str
    stderr: str
    success: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "stdout": self.stdout,
            "stderr": self.stderr,
            "success": self.success,
        }


@dataclass
class SessionInfo:
    session_id: str
    provider_name: str
    model: str
    created_at: str
    message_count: int = 0
