import { useEffect, useRef, useState, memo } from "react";
import { MessageBubble } from "./MessageBubble";
import { CodeExecutionCard } from "./CodeExecutionCard";
import { ThinkingBlock } from "./ThinkingBlock";
import { QuickActions } from "./QuickActions";
import type { ChatMessage } from "../api/types";

interface ThinkingState {
  content: string;
  isActive: boolean;
}

interface Props {
  messages: ChatMessage[];
  isStreaming: boolean;
  streamStartTime: number | null;
  thinking: ThinkingState;
  onQuickAction?: (prompt: string) => void;
}

function ElapsedTimer({ startTime }: { startTime: number }) {
  const [elapsed, setElapsed] = useState(0);

  useEffect(() => {
    const interval = setInterval(() => {
      setElapsed((Date.now() - startTime) / 1000);
    }, 100);
    return () => clearInterval(interval);
  }, [startTime]);

  const display =
    elapsed < 10
      ? `${elapsed.toFixed(1)}s`
      : `${Math.floor(elapsed / 60)}m ${Math.floor(elapsed % 60)}s`;

  return <span className="tabular-nums">{display}</span>;
}

function ChatPanelInner({ messages, isStreaming, streamStartTime, thinking, onQuickAction }: Props) {
  const bottomRef = useRef<HTMLDivElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;
    const isNearBottom =
      container.scrollHeight - container.scrollTop - container.clientHeight < 150;
    if (isNearBottom) {
      bottomRef.current?.scrollIntoView({ behavior: "smooth" });
    }
  }, [messages]);

  if (messages.length === 0) {
    return (
      <div className="flex-1 overflow-y-auto p-4 flex flex-col items-center justify-center text-[var(--text-secondary)]">
        <div className="text-4xl mb-4">🧬</div>
        <h2 className="text-xl font-semibold text-[var(--accent-cyan)] mb-2">
          EPIC
        </h2>
        <p className="text-sm text-center max-w-md mb-2">
          Engineering Probes via Instructing a Chatbot. Design QDB probes, PCR
          primers, and perform sequence analysis.
        </p>
        {onQuickAction && <QuickActions onAction={onQuickAction} />}
      </div>
    );
  }

  return (
    <div className="flex-1 overflow-y-auto p-4" ref={containerRef}>
      {messages.map((msg, i) => {
        // Skip empty non-streaming messages
        if (!msg.content && !msg.codeExecution && msg.isComplete) return null;

        if (msg.codeExecution) {
          return <CodeExecutionCard key={`code-${i}`} execution={msg.codeExecution} />;
        }

        // Skip empty assistant messages that haven't started streaming yet
        if (!msg.content && msg.role === "assistant" && !isStreaming) return null;
        if (!msg.content && msg.role === "assistant" && i < messages.length - 1) return null;

        return <MessageBubble key={i} message={msg} />;
      })}

      {(thinking.content || thinking.isActive) && (
        <ThinkingBlock content={thinking.content} isActive={thinking.isActive} />
      )}

      {isStreaming && streamStartTime && (
        <div className="flex items-center gap-2 text-[var(--text-secondary)] text-sm ml-2 mb-2">
          <div className="animate-pulse">●</div>
          <span>Generating...</span>
          <ElapsedTimer startTime={streamStartTime} />
        </div>
      )}
      <div ref={bottomRef} />
    </div>
  );
}

export const ChatPanel = memo(ChatPanelInner);
