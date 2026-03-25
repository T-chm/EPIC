interface Props {
  output: string;
  success: boolean;
}

export function CodeOutput({ output, success }: Props) {
  const borderColor = success ? "var(--accent-green)" : "var(--accent-red)";
  const label = success ? "Output" : "Error";

  return (
    <div
      className="my-2 rounded-lg overflow-hidden"
      style={{ border: `1px solid ${borderColor}` }}
    >
      <div
        className="px-3 py-1 text-xs font-medium"
        style={{ backgroundColor: borderColor, color: "var(--bg-primary)" }}
      >
        {label}
      </div>
      <pre className="p-3 text-sm bg-[var(--bg-secondary)] text-[var(--text-primary)] overflow-x-auto whitespace-pre-wrap">
        {output}
      </pre>
    </div>
  );
}
