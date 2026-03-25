import type { SessionInfo } from "../api/types";

interface Props {
  sessions: SessionInfo[];
  activeSessionId: string | null;
  onSelect: (id: string) => void;
  onCreate: () => void;
  onDelete: (id: string) => void;
}

export function Sidebar({
  sessions,
  activeSessionId,
  onSelect,
  onCreate,
  onDelete,
}: Props) {
  return (
    <div className="w-64 bg-[var(--bg-secondary)] border-r border-[var(--border-color)] flex flex-col h-full">
      <div className="p-4 border-b border-[var(--border-color)]">
        <h1 className="text-lg font-bold text-[var(--accent-cyan)]">EPIC</h1>
        <p className="text-xs text-[var(--text-secondary)]">
          Nucleic Acid Assay Design
        </p>
      </div>

      <div className="p-3">
        <button
          onClick={onCreate}
          className="w-full py-2 px-3 rounded-lg bg-[var(--accent-cyan)] text-[var(--bg-primary)] font-medium text-sm hover:opacity-90 transition-opacity"
        >
          + New Chat
        </button>
      </div>

      <div className="flex-1 overflow-y-auto px-3">
        {sessions.map((s) => (
          <div
            key={s.session_id}
            onClick={() => onSelect(s.session_id)}
            className={`group flex items-center justify-between rounded-lg px-3 py-2 mb-1 cursor-pointer text-sm transition-colors ${
              s.session_id === activeSessionId
                ? "bg-[var(--bg-tertiary)] text-[var(--text-primary)]"
                : "text-[var(--text-secondary)] hover:bg-[var(--bg-tertiary)]"
            }`}
          >
            <div className="min-w-0">
              <div className="truncate font-medium">
                {s.provider_name}/{s.model}
              </div>
              <div className="text-xs opacity-60">
                {s.message_count} messages
              </div>
            </div>
            <button
              onClick={(e) => {
                e.stopPropagation();
                onDelete(s.session_id);
              }}
              className="opacity-0 group-hover:opacity-100 text-[var(--accent-red)] hover:text-red-400 text-xs ml-2"
            >
              ✕
            </button>
          </div>
        ))}
      </div>
    </div>
  );
}
