interface Props {
  onAction: (prompt: string) => void;
}

const ACTIONS = [
  {
    label: "Design QDB Assay",
    icon: "🧬",
    description: "Design quantum dot barcode probes for a target sequence",
    prompt: "Design a nucleic acid quantum dot barcode assay probes for this sequence: ",
  },
  {
    label: "Design PCR Primers",
    icon: "🔬",
    description: "Design PCR primers for a template sequence",
    prompt: "Design PCR primers for the following template sequence: ",
  },
  {
    label: "Align Sequences",
    icon: "📊",
    description: "Multiple sequence alignment using MAFFT",
    prompt: "Align the following sequences using MAFFT:\n>seq1\n\n>seq2\n",
  },
  {
    label: "Analyze Probe",
    icon: "🔍",
    description: "Check GC content, Tm, and properties of a probe",
    prompt: "Analyze the following probe sequence for GC content, melting temperature, and other properties: ",
  },
];

export function QuickActions({ onAction }: Props) {
  return (
    <div className="grid grid-cols-2 gap-3 max-w-lg mt-6">
      {ACTIONS.map((action) => (
        <button
          key={action.label}
          onClick={() => onAction(action.prompt)}
          className="flex flex-col items-start gap-1 p-3 rounded-xl bg-[var(--bg-secondary)] border border-[var(--border-color)] hover:border-[var(--accent-cyan)] hover:bg-[var(--bg-tertiary)] transition-colors text-left"
        >
          <div className="flex items-center gap-2">
            <span className="text-lg">{action.icon}</span>
            <span className="text-sm font-medium text-[var(--text-primary)]">
              {action.label}
            </span>
          </div>
          <span className="text-xs text-[var(--text-secondary)] leading-tight">
            {action.description}
          </span>
        </button>
      ))}
    </div>
  );
}
