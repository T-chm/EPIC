import { useState, useCallback, useEffect } from "react";
import { Sidebar } from "./components/Sidebar";
import { ChatPanel } from "./components/ChatPanel";
import { InputBar } from "./components/InputBar";
import { ProviderSelector } from "./components/ProviderSelector";
import { SettingsPopover } from "./components/SettingsPopover";
import { useProviders } from "./hooks/useProviders";
import { useSessions } from "./hooks/useSessions";
import { useChat } from "./hooks/useChat";
import { switchProvider } from "./api/client";
import type { SessionInfo } from "./api/types";

export default function App() {
  const { providers } = useProviders();
  const { sessions, create, remove, refresh } = useSessions();
  const {
    messages,
    isStreaming,
    error,
    connected,
    streamStartTime,
    thinking,
    sendMessage,
    connect,
    disconnect,
    clearMessages,
  } = useChat();

  const [activeSession, setActiveSession] = useState<SessionInfo | null>(null);
  const [inputPrefill, setInputPrefill] = useState("");

  const handleCreateSession = useCallback(async () => {
    // Use first available provider
    const available = providers.find((p) => p.available);
    if (!available) return;
    const session = await create(available.name, available.default_model || undefined);
    setActiveSession(session);
    clearMessages();
    connect(session.session_id);
  }, [providers, create, clearMessages, connect]);

  const handleSelectSession = useCallback(
    (id: string) => {
      const session = sessions.find((s) => s.session_id === id);
      if (!session) return;
      disconnect();
      setActiveSession(session);
      clearMessages();
      connect(id);
    },
    [sessions, disconnect, clearMessages, connect]
  );

  const handleDeleteSession = useCallback(
    async (id: string) => {
      await remove(id);
      if (activeSession?.session_id === id) {
        disconnect();
        setActiveSession(null);
        clearMessages();
      }
    },
    [remove, activeSession, disconnect, clearMessages]
  );

  const handleSwitchProvider = useCallback(
    async (provider: string, model: string) => {
      if (!activeSession) return;
      try {
        const updated = await switchProvider(
          activeSession.session_id,
          provider,
          model
        );
        setActiveSession(updated);
        refresh();
      } catch (e) {
        console.error("Failed to switch provider:", e);
      }
    },
    [activeSession, refresh]
  );

  // Auto-create a session on first load if providers are available and no sessions exist
  useEffect(() => {
    if (providers.length > 0 && sessions.length === 0 && !activeSession) {
      handleCreateSession();
    }
  }, [providers, sessions.length, activeSession, handleCreateSession]);

  return (
    <div className="flex h-screen">
      <Sidebar
        sessions={sessions}
        activeSessionId={activeSession?.session_id || null}
        onSelect={handleSelectSession}
        onCreate={handleCreateSession}
        onDelete={handleDeleteSession}
      />

      <div className="flex-1 flex flex-col">
        {/* Header */}
        <div className="flex items-center justify-between px-4 py-3 border-b border-[var(--border-color)] bg-[var(--bg-secondary)]">
          <div className="flex items-center gap-3">
            <span className="text-sm font-medium text-[var(--text-secondary)]">
              {activeSession
                ? `Session: ${activeSession.session_id}`
                : "No active session"}
            </span>
            {activeSession && (
              <span
                className={`ml-2 inline-block w-2 h-2 rounded-full ${
                  connected ? "bg-[var(--accent-green)]" : "bg-[var(--accent-red)]"
                }`}
                title={connected ? "Connected" : "Disconnected"}
              />
            )}
          </div>
          <div className="flex items-center gap-2">
            {activeSession && (
              <ProviderSelector
                providers={providers}
                currentProvider={activeSession.provider_name}
                currentModel={activeSession.model}
                onSwitch={handleSwitchProvider}
                disabled={isStreaming}
              />
            )}
            <SettingsPopover disabled={isStreaming} />
          </div>
        </div>

        {/* Error banner */}
        {error && (
          <div className="px-4 py-2 bg-red-900/50 text-[var(--accent-red)] text-sm">
            {error}
          </div>
        )}

        {/* Chat area */}
        <ChatPanel
          messages={messages}
          isStreaming={isStreaming}
          streamStartTime={streamStartTime}
          thinking={thinking}
          onQuickAction={(prompt) => setInputPrefill(prompt)}
        />

        {/* Input */}
        <InputBar
          onSend={(msg) => { sendMessage(msg); setInputPrefill(""); }}
          disabled={isStreaming || !activeSession}
          prefill={inputPrefill}
        />
      </div>
    </div>
  );
}
