import { memo } from "react";
import ReactMarkdown from "react-markdown";
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import { oneDark } from "react-syntax-highlighter/dist/esm/styles/prism";
import type { ChatMessage } from "../api/types";

interface Props {
  message: ChatMessage;
}

function MessageBubbleInner({ message }: Props) {
  const isUser = message.role === "user";
  const isComplete = message.isComplete ?? isUser; // User messages are always complete

  return (
    <div className={`flex ${isUser ? "justify-end" : "justify-start"} mb-4`}>
      <div
        className={`max-w-[80%] rounded-xl px-4 py-3 ${
          isUser
            ? "bg-[var(--accent-blue)] text-white"
            : "bg-[var(--bg-secondary)] text-[var(--text-primary)]"
        }`}
      >
        <div className="text-xs font-medium mb-1 opacity-70">
          {isUser ? "You" : "Assistant"}
        </div>
        <div className="prose prose-invert prose-sm max-w-none">
          {isComplete ? (
            /* Full markdown rendering — only after streaming completes */
            <ReactMarkdown
              components={{
                code({ className, children, ...props }) {
                  const match = /language-(\w+)/.exec(className || "");
                  const inline = !match;
                  return inline ? (
                    <code
                      className="bg-[var(--bg-tertiary)] px-1 py-0.5 rounded text-sm"
                      {...props}
                    >
                      {children}
                    </code>
                  ) : (
                    <SyntaxHighlighter
                      style={oneDark}
                      language={match[1]}
                      customStyle={{
                        fontSize: "0.8rem",
                        borderRadius: "0.5rem",
                      }}
                    >
                      {String(children).replace(/\n$/, "")}
                    </SyntaxHighlighter>
                  );
                },
              }}
            >
              {message.content}
            </ReactMarkdown>
          ) : (
            /* Raw text during streaming — no markdown parsing overhead */
            <pre className="whitespace-pre-wrap font-sans text-sm m-0 leading-relaxed">
              {message.content}
              <span className="animate-pulse">▊</span>
            </pre>
          )}
        </div>
      </div>
    </div>
  );
}

export const MessageBubble = memo(MessageBubbleInner);
