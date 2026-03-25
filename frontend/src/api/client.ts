import type { ProviderInfo, SessionInfo } from "./types";

const BASE = "/api";

async function fetchJSON<T>(url: string, init?: RequestInit): Promise<T> {
  const res = await fetch(`${BASE}${url}`, {
    headers: { "Content-Type": "application/json" },
    ...init,
  });
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`${res.status}: ${body}`);
  }
  if (res.status === 204) return undefined as T;
  return res.json();
}

export async function getProviders(): Promise<ProviderInfo[]> {
  return fetchJSON("/providers");
}

export async function getProviderModels(name: string): Promise<{ provider: string; models: string[] }> {
  return fetchJSON(`/providers/${name}/models`);
}

export async function createSession(provider: string, model?: string): Promise<SessionInfo> {
  return fetchJSON("/sessions", {
    method: "POST",
    body: JSON.stringify({ provider, model }),
  });
}

export async function listSessions(): Promise<SessionInfo[]> {
  return fetchJSON("/sessions");
}

export async function getSession(id: string): Promise<SessionInfo & { messages: { role: string; content: string }[] }> {
  return fetchJSON(`/sessions/${id}`);
}

export async function deleteSession(id: string): Promise<void> {
  return fetchJSON(`/sessions/${id}`, { method: "DELETE" });
}

export async function switchProvider(sessionId: string, provider: string, model?: string): Promise<SessionInfo> {
  return fetchJSON(`/sessions/${sessionId}/provider`, {
    method: "PATCH",
    body: JSON.stringify({ provider, model }),
  });
}

export async function executeCode(sessionId: string, code: string): Promise<{ code: string; stdout: string; stderr: string; success: boolean }> {
  return fetchJSON("/code/execute", {
    method: "POST",
    body: JSON.stringify({ session_id: sessionId, code }),
  });
}
