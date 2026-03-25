import { useState, useEffect, useRef } from "react";

interface Props {
  disabled: boolean;
}

export function SettingsPopover({ disabled }: Props) {
  const [open, setOpen] = useState(false);
  const ref = useRef<HTMLDivElement>(null);

  const [theme, setTheme] = useState(() => localStorage.getItem("epic-theme") || "dark");
  const [autoExecute, setAutoExecute] = useState(() => localStorage.getItem("epic-auto-execute") !== "false");

  // Apply theme on mount and change
  useEffect(() => {
    document.documentElement.setAttribute("data-theme", theme);
    localStorage.setItem("epic-theme", theme);
  }, [theme]);

  useEffect(() => {
    localStorage.setItem("epic-auto-execute", String(autoExecute));
  }, [autoExecute]);

  // Close on outside click
  useEffect(() => {
    if (!open) return;
    const handler = (e: MouseEvent) => {
      if (ref.current && !ref.current.contains(e.target as Node)) {
        setOpen(false);
      }
    };
    document.addEventListener("mousedown", handler);
    return () => document.removeEventListener("mousedown", handler);
  }, [open]);

  return (
    <div className="relative" ref={ref}>
      <button
        onClick={() => setOpen(!open)}
        disabled={disabled}
        className="p-1.5 rounded hover:bg-[var(--bg-tertiary)] text-[var(--text-secondary)] transition-colors disabled:opacity-50"
        title="Settings"
      >
        <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
          <circle cx="12" cy="12" r="3" />
          <path d="M19.4 15a1.65 1.65 0 00.33 1.82l.06.06a2 2 0 010 2.83 2 2 0 01-2.83 0l-.06-.06a1.65 1.65 0 00-1.82-.33 1.65 1.65 0 00-1 1.51V21a2 2 0 01-4 0v-.09A1.65 1.65 0 009 19.4a1.65 1.65 0 00-1.82.33l-.06.06a2 2 0 01-2.83 0 2 2 0 010-2.83l.06-.06A1.65 1.65 0 004.68 15a1.65 1.65 0 00-1.51-1H3a2 2 0 010-4h.09A1.65 1.65 0 004.6 9a1.65 1.65 0 00-.33-1.82l-.06-.06a2 2 0 010-2.83 2 2 0 012.83 0l.06.06A1.65 1.65 0 009 4.68a1.65 1.65 0 001-1.51V3a2 2 0 014 0v.09a1.65 1.65 0 001 1.51 1.65 1.65 0 001.82-.33l.06-.06a2 2 0 012.83 0 2 2 0 010 2.83l-.06.06A1.65 1.65 0 0019.4 9a1.65 1.65 0 001.51 1H21a2 2 0 010 4h-.09a1.65 1.65 0 00-1.51 1z" />
        </svg>
      </button>

      {open && (
        <div className="absolute right-0 top-full mt-2 w-64 bg-[var(--bg-secondary)] border border-[var(--border-color)] rounded-xl shadow-xl z-50 p-4">
          <h3 className="text-sm font-semibold text-[var(--text-primary)] mb-3">
            Settings
          </h3>

          {/* Theme toggle */}
          <div className="flex items-center justify-between mb-3">
            <span className="text-sm text-[var(--text-secondary)]">Theme</span>
            <button
              onClick={() => setTheme(theme === "dark" ? "light" : "dark")}
              className="px-3 py-1 rounded-lg bg-[var(--bg-tertiary)] text-xs text-[var(--text-primary)] hover:opacity-80 transition-opacity"
            >
              {theme === "dark" ? "Dark" : "Light"}
            </button>
          </div>

          {/* Auto-execute toggle */}
          <div className="flex items-center justify-between mb-3">
            <span className="text-sm text-[var(--text-secondary)]">Auto-execute code</span>
            <button
              onClick={() => setAutoExecute(!autoExecute)}
              className={`w-10 h-5 rounded-full relative transition-colors ${
                autoExecute ? "bg-[var(--accent-green)]" : "bg-[var(--bg-tertiary)]"
              }`}
            >
              <span
                className={`absolute top-0.5 w-4 h-4 rounded-full bg-white transition-transform ${
                  autoExecute ? "translate-x-5" : "translate-x-0.5"
                }`}
              />
            </button>
          </div>

          <div className="text-xs text-[var(--text-secondary)] mt-2 pt-2 border-t border-[var(--border-color)]">
            Provider and model can be changed in the header bar.
          </div>
        </div>
      )}
    </div>
  );
}
