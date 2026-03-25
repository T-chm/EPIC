export interface ProviderInfo {
  name: string;
  available: boolean;
  models: string[];
  default_model: string | null;
}

export interface SessionInfo {
  session_id: string;
  provider_name: string;
  model: string;
  created_at: string;
  message_count: number;
}

export interface CodeExecution {
  code: string;
  status: "detected" | "executing" | "done";
  result?: string;
  success?: boolean;
}

export interface ChatMessage {
  role: "user" | "assistant";
  content: string;
  isComplete?: boolean;
  codeExecution?: CodeExecution;
}

export interface CodeBlock {
  code: string;
  result?: string;
  success?: boolean;
}

export interface WSTokenEvent {
  type: "token";
  content: string;
}

export interface WSCodeDetectedEvent {
  type: "code_detected";
  code: string;
}

export interface WSCodeExecutingEvent {
  type: "code_executing";
}

export interface WSCodeResultEvent {
  type: "code_result";
  stdout: string;
  success: boolean;
}

export interface WSDoneEvent {
  type: "done";
  full_response: string;
}

export interface WSErrorEvent {
  type: "error";
  message: string;
}

export interface WSThinkingStartEvent {
  type: "thinking_start";
}

export interface WSThinkingTokenEvent {
  type: "thinking_token";
  content: string;
}

export interface WSThinkingEndEvent {
  type: "thinking_end";
}

export type WSEvent =
  | WSTokenEvent
  | WSCodeDetectedEvent
  | WSCodeExecutingEvent
  | WSCodeResultEvent
  | WSDoneEvent
  | WSErrorEvent
  | WSThinkingStartEvent
  | WSThinkingTokenEvent
  | WSThinkingEndEvent;
