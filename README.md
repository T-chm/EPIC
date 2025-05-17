# EPIC: Engineering Probes via Instructing a Chatbot

EPIC is an LLM based chatbot designed specifically for nucleic acid diagnostic assay development. Currently, this platform supports the design, optimization, and troubleshooting of the quantum dot barcode (QDB) and PCR assays.

## Features

- **QDB Probe Design**: Design and optimize capture and reporter probes for QDB assays
- **PCR Primer Design**: Design optimized primers for PCR assays using Primer3
- **Multiple Sequence Alignment**: Built-in MAFFT integration for sequence alignment and target region identification
- **Thermodynamic Analysis**: Calculates optimal GC content, melting temperature, and other parameters
- **Expert Guidance**: Assay development advice following established design principles for both QDB and PCR

## Assay Design Processes

### QDB Assay Design

EPIC follows a systematic approach to QDB assay design:

1. Target sequence analysis and verification
2. Splitting the target sequence into equal-length halves 
3. Generating reverse complements for capture and reporter probes
4. Optimizing GC content (40-60%) and melting temperature (55-72°C)
5. Adding appropriate chemical modifications (/5AmMC6/ and /3Cy5Sp/)

### PCR Primer Design

For PCR primer design, EPIC follows these guidelines:

1. Target sequence analysis and amplicon region selection
2. Designing primers with optimal length (18-30 nucleotides)
3. Optimizing melting temperature (50-65°C) with less than 5°C difference between primers
4. Ensuring appropriate GC content (40-60%) with a G or C at the 3' end
5. Avoiding secondary structures, primer-dimers, and runs of 4+ identical nucleotides

## Requirements

- Python 3.9+
- OpenAI API key (GPT-4o model access required)
- MAFFT installed on your system (for sequence alignment features)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/proteinuniverse/epic-chatbot.git
cd epic-chatbot
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

3. Install MAFFT (if not already installed):
   - **macOS**: `brew install mafft`
   - **Ubuntu/Debian**: `sudo apt-get install mafft`
   - **Windows**: Download from the [MAFFT website](https://mafft.cbrc.jp/alignment/software/)

4. Create a `.env` file with your OpenAI API key:
```
OPENAI_API_KEY=your_api_key_here
```

## Usage

Run the EPIC chatbot:
```bash
python epic_chat.py
```

This will start a Gradio web interface accessible at http://localhost:7860.

### Example Queries

The repository includes sample queries in the `data/sample_user_prompts.txt` file. Here are some examples:

#### QDB Assay Queries
- "Design QDB probes for detecting this sequence: [paste your sequence]"
- "Analyze these variant sequences for conserved regions"
- "How do I optimize my QDB assay for higher sensitivity?"
- "Design a QDB assay for detecting HIV variants"

#### PCR Assay Queries
- "Design PCR primers for amplifying this region: [paste your sequence]"
- "Create primers for a 200bp amplicon from this sequence"
- "What are the best primer design parameters for a multiplex PCR?"
- "Optimize these primers for my SARS-CoV-2 detection assay"

### Data Folder

The `data` folder contains:
- Sample user prompts for reference in `sample_user_prompts.txt`
- Reference materials for QDB assay design
- Example sequences for testing and demonstration

### Python Code Execution

You can execute Python code in two ways:
1. Prefix your message with `!code` followed by your Python code
2. EPIC will automatically execute any Python code blocks it generates in its responses

## How It Works

EPIC leverages the OpenAI GPT-4o model with specialized training in molecular biology and nucleic acid diagnostics. The system includes:

1. Specialized prompt engineering for both QDB and PCR assay design
2. Integration with BioPython for sequence analysis and thermodynamic calculations
3. MAFFT integration for multiple sequence alignment of variant sequences
4. Interactive code execution for real-time testing and visualization
5. Built-in primer design using primer3-py for optimized PCR primer generation
6. Step-by-step assay optimization guidance based on best practices

## License

This project is licensed under the Educational Community License Version 2.0 (ECL-2.0) - see the [LICENSE](LICENSE) file for details.

## Acknowledgements

- [Gradio](https://gradio.app/) for the web interface
- [OpenAI](https://openai.com/) for the GPT-4o model
- [BioPython](https://biopython.org/) for biological sequence analysis
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) for multiple sequence alignment
- [Primer3-py](https://libnano.github.io/primer3-py/) for PCR primer design
