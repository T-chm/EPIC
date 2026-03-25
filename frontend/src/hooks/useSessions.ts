import { useState, useEffect, useCallback } from "react";
import {
  listSessions,
  createSession,
  deleteSession as apiDeleteSession,
} from "../api/client";
import type { SessionInfo } from "../api/types";

export function useSessions() {
  const [sessions, setSessions] = useState<SessionInfo[]>([]);
  const [loading, setLoading] = useState(true);

  const refresh = useCallback(async () => {
    try {
      const data = await listSessions();
      setSessions(data);
    } catch {
      // silently fail
    } finally {
      setLoading(false);
    }
  }, []);

  useEffect(() => {
    refresh();
  }, [refresh]);

  const create = useCallback(
    async (provider: string, model?: string): Promise<SessionInfo> => {
      const session = await createSession(provider, model);
      setSessions((prev) => [...prev, session]);
      return session;
    },
    []
  );

  const remove = useCallback(
    async (id: string) => {
      await apiDeleteSession(id);
      setSessions((prev) => prev.filter((s) => s.session_id !== id));
    },
    []
  );

  return { sessions, loading, create, remove, refresh };
}
