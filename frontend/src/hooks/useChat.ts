import { useState, useRef, useCallback, useEffect } from "react";
import type { ChatMessage, CodeBlock, WSEvent } from "../api/types";

interface ThinkingState {
  content: string;
  isActive: boolean;
}

interface UseChatReturn {
  messages: ChatMessage[];
  codeBlocks: CodeBlock[];
  isStreaming: boolean;
  error: string | null;
  connected: boolean;
  streamStartTime: number | null;
  thinking: ThinkingState;
  sendMessage: (content: string) => void;
  connect: (sessionId: string) => void;
  disconnect: () => void;
  clearMessages: () => void;
}

export function useChat(): UseChatReturn {
  const [messages, setMessages] = useState<ChatMessage[]>([]);
  const [codeBlocks, setCodeBlocks] = useState<CodeBlock[]>([]);
  const [isStreaming, setIsStreaming] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [connected, setConnected] = useState(false);
  const [streamStartTime, setStreamStartTime] = useState<number | null>(null);
  const [thinking, setThinking] = useState<ThinkingState>({ content: "", isActive: false });
  const wsRef = useRef<WebSocket | null>(null);

  // Token buffering via requestAnimationFrame
  const tokenBufferRef = useRef("");
  const rafRef = useRef<number | null>(null);

  const flushTokenBuffer = useCallback(() => {
    rafRef.current = null;
    const buffered = tokenBufferRef.current;
    if (!buffered) return;
    tokenBufferRef.current = "";

    setMessages((prev) => {
      const updated = [...prev];
      const last = updated[updated.length - 1];
      if (last && last.role === "assistant" && !last.codeExecution) {
        updated[updated.length - 1] = {
          ...last,
          content: last.content + buffered,
        };
      }
      return updated;
    });
  }, []);

  const scheduleFlush = useCallback(() => {
    if (rafRef.current === null) {
      rafRef.current = requestAnimationFrame(flushTokenBuffer);
    }
  }, [flushTokenBuffer]);

  useEffect(() => {
    return () => {
      if (rafRef.current !== null) cancelAnimationFrame(rafRef.current);
    };
  }, []);

  const connect = useCallback((sessionId: string) => {
    if (wsRef.current) {
      wsRef.current.onclose = null;
      wsRef.current.close();
      wsRef.current = null;
    }

    const protocol = window.location.protocol === "https:" ? "wss:" : "ws:";
    const ws = new WebSocket(
      `${protocol}//${window.location.host}/ws/chat/${sessionId}`
    );

    ws.onopen = () => {
      setError(null);
      setConnected(true);
    };

    ws.onmessage = (event) => {
      let data: WSEvent;
      try {
        data = JSON.parse(event.data);
      } catch {
        return;
      }

      switch (data.type) {
        case "token":
          tokenBufferRef.current += data.content;
          scheduleFlush();
          break;

        case "thinking_start":
          setThinking({ content: "", isActive: true });
          break;

        case "thinking_token":
          setThinking((prev) => ({ ...prev, content: prev.content + data.content }));
          break;

        case "thinking_end":
          setThinking((prev) => ({ ...prev, isActive: false }));
          break;

        case "code_detected":
          // Flush any pending tokens first
          if (tokenBufferRef.current) {
            if (rafRef.current !== null) {
              cancelAnimationFrame(rafRef.current);
              rafRef.current = null;
            }
            const remaining = tokenBufferRef.current;
            tokenBufferRef.current = "";
            setMessages((prev) => {
              const updated = [...prev];
              const last = updated[updated.length - 1];
              if (last && last.role === "assistant" && !last.codeExecution) {
                updated[updated.length - 1] = { ...last, content: last.content + remaining };
              }
              return updated;
            });
          }
          // Mark current assistant message as complete before code card
          setMessages((prev) => {
            const updated = [...prev];
            const last = updated[updated.length - 1];
            if (last && last.role === "assistant" && !last.codeExecution && last.content) {
              updated[updated.length - 1] = { ...last, isComplete: true };
            }
            return [
              ...updated,
              {
                role: "assistant",
                content: "",
                codeExecution: { code: data.code, status: "detected" },
              },
            ];
          });
          setCodeBlocks((prev) => [...prev, { code: data.code }]);
          break;

        case "code_executing":
          setMessages((prev) => {
            const updated = [...prev];
            for (let i = updated.length - 1; i >= 0; i--) {
              if (updated[i].codeExecution && updated[i].codeExecution!.status === "detected") {
                updated[i] = {
                  ...updated[i],
                  codeExecution: { ...updated[i].codeExecution!, status: "executing" },
                };
                break;
              }
            }
            return updated;
          });
          break;

        case "code_result":
          setMessages((prev) => {
            const updated = [...prev];
            for (let i = updated.length - 1; i >= 0; i--) {
              if (updated[i].codeExecution && updated[i].codeExecution!.status !== "done") {
                updated[i] = {
                  ...updated[i],
                  codeExecution: {
                    ...updated[i].codeExecution!,
                    status: "done",
                    result: data.stdout,
                    success: data.success,
                  },
                };
                break;
              }
            }
            // Add an empty assistant message for the interpretation turn
            return [...updated, { role: "assistant" as const, content: "" }];
          });
          setCodeBlocks((prev) => {
            const updated = [...prev];
            const last = updated[updated.length - 1];
            if (last && !last.result) {
              updated[updated.length - 1] = { ...last, result: data.stdout, success: data.success };
            }
            return updated;
          });
          break;

        case "done":
          // Flush remaining tokens
          if (tokenBufferRef.current) {
            if (rafRef.current !== null) {
              cancelAnimationFrame(rafRef.current);
              rafRef.current = null;
            }
            const remaining = tokenBufferRef.current;
            tokenBufferRef.current = "";
            setMessages((prev) => {
              const updated = [...prev];
              const last = updated[updated.length - 1];
              if (last && last.role === "assistant" && !last.codeExecution) {
                updated[updated.length - 1] = { ...last, content: last.content + remaining };
              }
              return updated;
            });
          }
          // Mark all assistant messages as complete
          setMessages((prev) =>
            prev.map((m) =>
              m.role === "assistant" && !m.isComplete && !m.codeExecution
                ? { ...m, isComplete: true }
                : m
            )
          );
          setIsStreaming(false);
          setStreamStartTime(null);
          break;

        case "error":
          setError(data.message);
          setIsStreaming(false);
          setStreamStartTime(null);
          break;
      }
    };

    ws.onerror = () => {
      setError("WebSocket connection error");
      setIsStreaming(false);
      setStreamStartTime(null);
    };

    ws.onclose = (ev) => {
      setConnected(false);
      if (ev.code !== 1000 && wsRef.current === ws) {
        setError("Connection lost. Click 'New Chat' to reconnect.");
      }
      setIsStreaming(false);
    };

    wsRef.current = ws;
  }, [scheduleFlush]);

  const disconnect = useCallback(() => {
    if (wsRef.current) {
      wsRef.current.onclose = null;
      wsRef.current.close();
      wsRef.current = null;
    }
    setConnected(false);
  }, []);

  const sendMessage = useCallback((content: string) => {
    if (!wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) {
      setError("Not connected. Create or select a session first.");
      return;
    }

    setMessages((prev) => [
      ...prev,
      { role: "user", content, isComplete: true },
      { role: "assistant", content: "" },
    ]);
    setIsStreaming(true);
    setStreamStartTime(Date.now());
    setError(null);
    setThinking({ content: "", isActive: false });
    tokenBufferRef.current = "";

    wsRef.current.send(JSON.stringify({ type: "message", content }));
  }, []);

  const clearMessages = useCallback(() => {
    setMessages([]);
    setCodeBlocks([]);
    setError(null);
    tokenBufferRef.current = "";
    setStreamStartTime(null);
  }, []);

  return {
    messages,
    codeBlocks,
    isStreaming,
    error,
    connected,
    streamStartTime,
    thinking,
    sendMessage,
    connect,
    disconnect,
    clearMessages,
  };
}
