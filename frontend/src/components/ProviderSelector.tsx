import { useState, useEffect } from "react";
import type { ProviderInfo } from "../api/types";

interface Props {
  providers: ProviderInfo[];
  currentProvider: string;
  currentModel: string;
  onSwitch: (provider: string, model: string) => void;
  disabled: boolean;
}

export function ProviderSelector({
  providers,
  currentProvider,
  currentModel,
  onSwitch,
  disabled,
}: Props) {
  const [selectedProvider, setSelectedProvider] = useState(currentProvider);
  const [selectedModel, setSelectedModel] = useState(currentModel);

  const availableProviders = providers.filter((p) => p.available);
  const currentProviderInfo = providers.find((p) => p.name === selectedProvider);
  const models = currentProviderInfo?.models || [];

  useEffect(() => {
    setSelectedProvider(currentProvider);
    setSelectedModel(currentModel);
  }, [currentProvider, currentModel]);

  const handleProviderChange = (name: string) => {
    setSelectedProvider(name);
    const provider = providers.find((p) => p.name === name);
    const defaultModel = provider?.default_model || provider?.models[0] || "";
    setSelectedModel(defaultModel);
    onSwitch(name, defaultModel);
  };

  const handleModelChange = (model: string) => {
    setSelectedModel(model);
    onSwitch(selectedProvider, model);
  };

  return (
    <div className="flex items-center gap-3">
      <select
        value={selectedProvider}
        onChange={(e) => handleProviderChange(e.target.value)}
        disabled={disabled}
        className="bg-[var(--bg-tertiary)] text-[var(--text-primary)] rounded px-3 py-1.5 text-sm outline-none border border-[var(--border-color)] disabled:opacity-50"
      >
        {availableProviders.map((p) => (
          <option key={p.name} value={p.name}>
            {p.name}
          </option>
        ))}
      </select>

      <select
        value={selectedModel}
        onChange={(e) => handleModelChange(e.target.value)}
        disabled={disabled}
        className="bg-[var(--bg-tertiary)] text-[var(--text-primary)] rounded px-3 py-1.5 text-sm outline-none border border-[var(--border-color)] disabled:opacity-50"
      >
        {models.map((m) => (
          <option key={m} value={m}>
            {m}
          </option>
        ))}
      </select>
    </div>
  );
}
