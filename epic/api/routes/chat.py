"""WebSocket chat endpoint for streaming LLM responses."""

from __future__ import annotations

import json
import logging
import traceback

from fastapi import APIRouter, WebSocket, WebSocketDisconnect
from starlette.websockets import WebSocketState

from epic.api.dependencies import get_engine

logger = logging.getLogger(__name__)

router = APIRouter()


@router.websocket("/ws/chat/{session_id}")
async def chat_websocket(websocket: WebSocket, session_id: str) -> None:
    """WebSocket endpoint for streaming chat."""
    await websocket.accept()
    logger.info("WebSocket connected for session %s", session_id)

    engine = get_engine(session_id)
    if not engine:
        await websocket.send_json({"type": "error", "message": "Session not found"})
        await websocket.close(code=4004, reason="Session not found")
        return

    try:
        while True:
            raw = await websocket.receive_text()
            try:
                data = json.loads(raw)
            except json.JSONDecodeError:
                await _safe_send(websocket, {"type": "error", "message": "Invalid JSON"})
                continue

            msg_type = data.get("type")
            if msg_type != "message":
                await _safe_send(websocket, {"type": "error", "message": f"Unknown type: {msg_type}"})
                continue

            content = data.get("content", "").strip()
            if not content:
                await _safe_send(websocket, {"type": "error", "message": "Empty message"})
                continue

            # Stream response events to the client
            try:
                async for event in engine.astream_response(content):
                    if event.event_type == "token":
                        await _safe_send(websocket, {"type": "token", "content": event.data})
                    elif event.event_type == "thinking_start":
                        await _safe_send(websocket, {"type": "thinking_start"})
                    elif event.event_type == "thinking_token":
                        await _safe_send(websocket, {"type": "thinking_token", "content": event.data})
                    elif event.event_type == "thinking_end":
                        await _safe_send(websocket, {"type": "thinking_end"})
                    elif event.event_type == "code_detected":
                        await _safe_send(websocket, {"type": "code_detected", "code": event.data})
                    elif event.event_type == "code_executing":
                        await _safe_send(websocket, {"type": "code_executing"})
                    elif event.event_type == "code_result":
                        is_error = event.data.startswith("Error:")
                        await _safe_send(websocket, {
                            "type": "code_result",
                            "stdout": event.data,
                            "success": not is_error,
                        })
                    elif event.event_type == "done":
                        await _safe_send(websocket, {
                            "type": "done",
                            "full_response": event.data,
                        })
                    elif event.event_type == "error":
                        await _safe_send(websocket, {"type": "error", "message": event.data})
            except Exception as e:
                logger.error("Stream error for session %s: %s\n%s", session_id, e, traceback.format_exc())
                await _safe_send(websocket, {"type": "error", "message": f"Stream error: {e}"})
                # Don't break — keep the connection alive for the next message

    except WebSocketDisconnect:
        logger.info("WebSocket disconnected for session %s", session_id)
    except Exception as e:
        logger.error("WebSocket fatal error for session %s: %s\n%s", session_id, e, traceback.format_exc())
        await _safe_send(websocket, {"type": "error", "message": str(e)})


async def _safe_send(websocket: WebSocket, data: dict) -> None:
    """Send JSON over WebSocket, silently ignoring if connection is closed."""
    try:
        if websocket.client_state == WebSocketState.CONNECTED:
            await websocket.send_json(data)
    except Exception:
        pass
