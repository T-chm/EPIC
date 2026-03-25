import { memo, useState } from "react";
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import { oneDark } from "react-syntax-highlighter/dist/esm/styles/prism";
import type { CodeExecution } from "../api/types";

interface Props {
  execution: CodeExecution;
}

function CodeExecutionCardInner({ execution }: Props) {
  const [codeExpanded, setCodeExpanded] = useState(execution.status !== "done");
  const isDone = execution.status === "done";
  const isError = isDone && !execution.success;

  const statusIcon =
    execution.status === "detected"
      ? "◆"
      : execution.status === "executing"
        ? "⟳"
        : execution.success
          ? "✓"
          : "✗";

  const statusColor =
    execution.status === "executing"
      ? "var(--accent-cyan)"
      : isDone && execution.success
        ? "var(--accent-green)"
        : isDone
          ? "var(--accent-red)"
          : "var(--text-secondary)";

  const statusLabel =
    execution.status === "detected"
      ? "Code detected"
      : execution.status === "executing"
        ? "Executing..."
        : execution.success
          ? "Execution complete"
          : "Execution failed";

  return (
    <div className="mb-4 flex justify-start">
      <div className="max-w-[85%] w-full rounded-xl overflow-hidden border border-[var(--border-color)] bg-[var(--bg-secondary)]">
        {/* Header */}
        <div
          className="flex items-center justify-between px-4 py-2 cursor-pointer select-none"
          style={{ borderBottom: `1px solid var(--border-color)` }}
          onClick={() => setCodeExpanded(!codeExpanded)}
        >
          <div className="flex items-center gap-2">
            <span
              className={`text-sm font-bold ${execution.status === "executing" ? "animate-spin" : ""}`}
              style={{ color: statusColor }}
            >
              {statusIcon}
            </span>
            <span className="text-xs font-medium" style={{ color: statusColor }}>
              {statusLabel}
            </span>
            <span className="text-xs text-[var(--text-secondary)]">Python</span>
          </div>
          <span className="text-xs text-[var(--text-secondary)]">
            {codeExpanded ? "▼" : "▶"} {execution.code.split("\n").length} lines
          </span>
        </div>

        {/* Code block (collapsible) */}
        {codeExpanded && (
          <SyntaxHighlighter
            language="python"
            style={oneDark}
            customStyle={{
              margin: 0,
              fontSize: "0.75rem",
              maxHeight: "300px",
              background: "var(--bg-primary)",
            }}
            showLineNumbers
          >
            {execution.code}
          </SyntaxHighlighter>
        )}

        {/* Output */}
        {isDone && execution.result && (
          <div
            className="border-t"
            style={{ borderColor: statusColor }}
          >
            <div
              className="px-3 py-1 text-xs font-medium"
              style={{ backgroundColor: statusColor, color: "var(--bg-primary)" }}
            >
              {isError ? "Error" : "Output"}
            </div>
            <pre className="px-3 py-2 text-xs overflow-x-auto whitespace-pre-wrap m-0 max-h-[200px] overflow-y-auto">
              {execution.result}
            </pre>
          </div>
        )}

        {/* Loading bar for executing state */}
        {execution.status === "executing" && (
          <div className="h-0.5 bg-[var(--bg-tertiary)] overflow-hidden">
            <div
              className="h-full w-1/3 rounded"
              style={{
                background: "var(--accent-cyan)",
                animation: "slideRight 1.5s ease-in-out infinite",
              }}
            />
          </div>
        )}
      </div>
    </div>
  );
}

export const CodeExecutionCard = memo(CodeExecutionCardInner);
