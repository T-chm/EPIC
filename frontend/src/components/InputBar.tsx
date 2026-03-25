import { useState, useCallback, useEffect, type KeyboardEvent } from "react";

interface Props {
  onSend: (message: string) => void;
  disabled: boolean;
  prefill?: string;
}

export function InputBar({ onSend, disabled, prefill }: Props) {
  const [input, setInput] = useState("");

  // When prefill changes (from quick actions), populate the input
  useEffect(() => {
    if (prefill) {
      setInput(prefill);
    }
  }, [prefill]);

  const handleSend = useCallback(() => {
    const trimmed = input.trim();
    if (!trimmed || disabled) return;
    onSend(trimmed);
    setInput("");
  }, [input, disabled, onSend]);

  const handleKeyDown = (e: KeyboardEvent<HTMLTextAreaElement>) => {
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      handleSend();
    }
  };

  return (
    <div className="border-t border-[var(--border-color)] p-4 bg-[var(--bg-secondary)]">
      <div className="flex gap-3 max-w-4xl mx-auto">
        <textarea
          value={input}
          onChange={(e) => setInput(e.target.value)}
          onKeyDown={handleKeyDown}
          placeholder="Ask about QDB assay design, PCR primers, or type !code for direct execution..."
          disabled={disabled}
          rows={2}
          className="flex-1 bg-[var(--bg-tertiary)] text-[var(--text-primary)] rounded-lg px-4 py-3 resize-none outline-none placeholder:text-[var(--text-secondary)] focus:ring-2 focus:ring-[var(--accent-cyan)] disabled:opacity-50"
        />
        <button
          onClick={handleSend}
          disabled={disabled || !input.trim()}
          className="px-6 py-3 bg-[var(--accent-cyan)] text-[var(--bg-primary)] font-semibold rounded-lg hover:opacity-90 transition-opacity disabled:opacity-30"
        >
          Send
        </button>
      </div>
    </div>
  );
}
