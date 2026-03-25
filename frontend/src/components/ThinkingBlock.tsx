import { memo, useState } from "react";

interface Props {
  content: string;
  isActive: boolean;
}

function ThinkingBlockInner({ content, isActive }: Props) {
  const [expanded, setExpanded] = useState(false);

  if (!content && !isActive) return null;

  return (
    <div className="mb-4 flex justify-start">
      <div
        className="max-w-[80%] rounded-xl overflow-hidden border border-[var(--border-color)] bg-[var(--bg-secondary)]"
        style={{ opacity: 0.7 }}
      >
        <div
          className="flex items-center gap-2 px-3 py-1.5 cursor-pointer select-none text-xs"
          onClick={() => setExpanded(!expanded)}
        >
          <span className={isActive ? "animate-pulse" : ""}>
            {isActive ? "💭" : "💭"}
          </span>
          <span className="text-[var(--text-secondary)] font-medium">
            {isActive ? "Thinking..." : "Thought"}
          </span>
          <span className="text-[var(--text-secondary)]">
            {expanded ? "▼" : "▶"}
          </span>
        </div>
        {expanded && (
          <div className="px-3 pb-2 max-h-[200px] overflow-y-auto">
            <pre className="text-xs text-[var(--text-secondary)] italic whitespace-pre-wrap font-sans m-0 leading-relaxed">
              {content}
            </pre>
          </div>
        )}
      </div>
    </div>
  );
}

export const ThinkingBlock = memo(ThinkingBlockInner);
