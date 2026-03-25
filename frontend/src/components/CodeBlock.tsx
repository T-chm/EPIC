import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import { oneDark } from "react-syntax-highlighter/dist/esm/styles/prism";

interface Props {
  code: string;
  onRun?: (code: string) => void;
}

export function CodeBlock({ code, onRun }: Props) {
  return (
    <div className="my-2 rounded-lg overflow-hidden border border-[var(--border-color)]">
      <div className="flex items-center justify-between bg-[var(--bg-tertiary)] px-3 py-1.5">
        <span className="text-xs text-[var(--text-secondary)]">Python</span>
        {onRun && (
          <button
            onClick={() => onRun(code)}
            className="text-xs px-2 py-0.5 rounded bg-[var(--accent-green)] text-[var(--bg-primary)] hover:opacity-80 transition-opacity"
          >
            Run
          </button>
        )}
      </div>
      <SyntaxHighlighter
        language="python"
        style={oneDark}
        customStyle={{
          margin: 0,
          fontSize: "0.8rem",
          background: "var(--bg-secondary)",
        }}
      >
        {code}
      </SyntaxHighlighter>
    </div>
  );
}
