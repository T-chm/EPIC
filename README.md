# EPIC: Engineering Probes via Instructing a Chatbot

EPIC is a specialized AI assistant designed specifically for quantum dot barcode (QDB) assay development. This platform helps researchers design, optimize, and troubleshoot QDB assays for nucleic acid detection.

![EPIC Chatbot](assets/epic_logo.png)


## Features

- **QDB Probe Design**: Design and optimize capture and reporter probes for QDB assays
- **QDB Assay Development**: Expert guidance on QDB assay development following established design principles
- **Multiple Sequence Alignment**: Built-in MAFFT integration for sequence alignment and target region identification
- **Thermodynamic Analysis**: Calculates optimal GC content, melting temperature, and other parameters for QDB probes

## QDB Assay Design Process

EPIC follows a systematic approach to QDB assay design:

1. Target sequence analysis and verification
2. Splitting the target sequence into equal-length halves 
3. Generating reverse complements for capture and reporter probes
4. Optimizing GC content (40-60%) and melting temperature (55-72°C)
5. Adding appropriate chemical modifications (/5AmMC6/ and /3Cy5Sp/)

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

- "Design QDB probes for detecting this sequence: [paste your sequence]"
- "Analyze these variant sequences for conserved regions"
- "How do I optimize my QDB assay for higher sensitivity?"
- "Design a QDB assay for detecting HIV variants"

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

EPIC leverages the OpenAI GPT-4o model with specialized training in biochemistry and QDB assay design. The system includes:

1. A specialized prompt engineering for QDB assay design
2. Integration with BioPython for sequence analysis and thermodynamic calculations
3. MAFFT integration for multiple sequence alignment
4. Interactive code execution for real-time testing
5. Built-in primer design using primer3-py

## License

This project is licensed under the Educational Community License Version 2.0 (ECL-2.0) - see the [LICENSE](LICENSE) file for details.

## Acknowledgements

- [Gradio](https://gradio.app/) for the web interface
- [OpenAI](https://openai.com/) for the GPT-4o model
- [BioPython](https://biopython.org/) for biological sequence analysis
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) for multiple sequence alignment
