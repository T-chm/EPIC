"""Code execution endpoint."""

from __future__ import annotations

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from epic.api.dependencies import get_engine

router = APIRouter()


class ExecuteCodeRequest(BaseModel):
    code: str
    session_id: str


@router.post("/execute")
async def execute_code(req: ExecuteCodeRequest) -> dict:
    """Execute Python code within a session's interpreter."""
    engine = get_engine(req.session_id)
    if not engine:
        raise HTTPException(status_code=404, detail="Session not found")

    result = engine.execute_code(req.code)
    return result.to_dict()
