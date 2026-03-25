# EPIC: Engineering Probes via Instructing a Chatbot

EPIC is an AI-powered chatbot for nucleic acid diagnostic assay development. It supports the design, optimization, and troubleshooting of quantum dot barcode (QDB) and PCR assays using multiple LLM providers.

## Features

- **Multi-Provider LLM Support**: OpenAI, Anthropic Claude, and Ollama (local models) — each using native SDKs
- **QDB Probe Design**: Design and optimize capture and reporter probes for QDB assays
- **PCR Primer Design**: Optimized primers for PCR assays using Primer3
- **Multiple Sequence Alignment**: Built-in MAFFT integration
- **Automatic Code Execution**: LLM-generated Python code is auto-executed with BioPython and Primer3
- **Web UI**: FastAPI backend + React frontend with streaming chat
- **CLI Tool**: Interactive chat mode and batch FASTA processing

## Requirements

- Python 3.11+
- Node.js 18+ (for frontend development)
- MAFFT (optional, for sequence alignment)

## Installation

```bash
git clone https://github.com/T-chm/EPIC.git
cd EPIC
```

Install with your preferred provider(s):

```bash
# Ollama only (local models, no API key needed)
pip install -e ".[ollama,cli]"

# All providers
pip install -e ".[all-providers,cli]"

# Development
pip install -e ".[all-providers,cli,dev]"
```

Install MAFFT (optional):
- **macOS**: `brew install mafft`
- **Ubuntu/Debian**: `sudo apt-get install mafft`

## Configuration

Copy the example environment file and configure:

```bash
cp .env.example .env
```

Key settings:
```
DEFAULT_PROVIDER=ollama          # openai | anthropic | ollama
OLLAMA_HOST=http://localhost:11434
OLLAMA_DEFAULT_MODEL=qwen3.5:4b
OPENAI_API_KEY=sk-...           # if using OpenAI
ANTHROPIC_API_KEY=sk-ant-...    # if using Anthropic
```

## Usage

### CLI — Interactive Chat

```bash
epic chat --provider ollama --model qwen3.5:4b
epic chat --provider anthropic --model claude-sonnet-4-20250514
epic chat --provider openai --model gpt-4o
```

Interactive commands:
- `/provider <name>` — switch provider
- `/model <name>` — switch model
- `/history` — show conversation
- `/clear` — clear history
- `!code <python>` — execute code directly
- `/exit` — quit

### CLI — Batch Processing

```bash
epic batch sequences.fasta --provider ollama --pipeline qdb --output ./results
epic batch sequences.fasta --pipeline pcr --output ./results
epic batch sequences.fasta --pipeline full --output ./results
```

### CLI — List Providers

```bash
epic providers
```

### Web UI

Start the API server:

```bash
epic serve --port 8000
```

For development with hot-reload frontend:

```bash
# Terminal 1: API server
epic serve

# Terminal 2: React dev server
cd frontend
npm install
npm run dev
```

For production, build the frontend and it will be served by FastAPI:

```bash
cd frontend && npm run build
epic serve
```

## API Endpoints

| Method | Path | Description |
|--------|------|-------------|
| `GET` | `/api/health` | Health check |
| `GET` | `/api/providers` | List available providers and models |
| `POST` | `/api/sessions` | Create a chat session |
| `GET` | `/api/sessions` | List sessions |
| `DELETE` | `/api/sessions/{id}` | Delete a session |
| `PATCH` | `/api/sessions/{id}/provider` | Switch provider mid-session |
| `POST` | `/api/code/execute` | Execute Python code |
| `WS` | `/ws/chat/{session_id}` | Streaming chat via WebSocket |

## Assay Design

### QDB Assay Design

1. Target sequence analysis and verification
2. Splitting into equal-length halves
3. Generating reverse complements for capture and reporter probes
4. Optimizing GC content (40-60%) and melting temperature (55-72°C)
5. Adding chemical modifications (/5AmMC6/ and /3Cy5Sp/)

### PCR Primer Design

1. Amplicon region selection
2. Primer length optimization (18-30 nucleotides)
3. Melting temperature optimization (50-65°C, within 5°C between primers)
4. GC content optimization (40-60%) with G/C at 3' end
5. Avoiding secondary structures and primer-dimers

## Project Structure

```
EPIC/
├── epic/                  # Python package
│   ├── config.py          # Settings (pydantic-settings)
│   ├── models.py          # Data models
│   ├── chat.py            # ChatEngine (shared core)
│   ├── prompts.py         # System prompt
│   ├── providers/         # LLM providers (OpenAI, Anthropic, Ollama)
│   ├── tools/             # Interpreter, MAFFT
│   ├── api/               # FastAPI backend
│   └── cli/               # Typer CLI
├── frontend/              # React + Vite + Tailwind
├── tests/                 # Test suite
├── scripts/               # Batch analysis scripts
└── data/                  # Sample data
```

## License

Educational Community License Version 2.0 (ECL-2.0) — see [LICENSE](LICENSE).

## Acknowledgements

- [BioPython](https://biopython.org/) for biological sequence analysis
- [Primer3-py](https://libnano.github.io/primer3-py/) for PCR primer design
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) for multiple sequence alignment
- [FastAPI](https://fastapi.tiangolo.com/) for the web API
- [Ollama](https://ollama.com/) for local model inference
