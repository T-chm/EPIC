#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EPIC Chatbot: Enhanced Platform for Innovative Chemistry

A specialized AI assistant designed specifically for quantum dot barcode (QDB) 
assay development. This platform helps researchers design, optimize, and troubleshoot 
QDB assays for nucleic acid detection.
"""

# Standard library imports
import os
import sys
import io
import re
import time
import code
import tempfile
import subprocess
from pathlib import Path

# Third-party imports
import gradio as gr
import openai
import primer3
from dotenv import load_dotenv
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqUtils

# Initialize environment
load_dotenv()

# Initialize OpenAI API with proper error handling
api_key = os.getenv("OPENAI_API_KEY", "")
if api_key:
    api_key = api_key.strip('"\'')
    openai.api_key = api_key

def chat_with_gpt4(history):
    """
    Function to communicate with the OpenAI GPT-4 API and stream responses.
    
    Args:
        history (list): A list of message dictionaries containing conversation history.
        
    Yields:
        str: Chunks of the response from the GPT-4 model.
    """
    try:
        print("Starting chat completion with streaming...")
        
        # Create streaming response from OpenAI
        response_stream = openai.chat.completions.create(
            model="gpt-4o",
            messages=history,
            temperature=0.2,
            stream=True
        )

        # Process each chunk of the streamed response
        for chunk in response_stream:
            if chunk.choices[0].delta and chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content
                
    except Exception as e:
        error_msg = str(e)
        print(f"Error in chat_with_gpt4: {error_msg}")
        yield f"Error: {error_msg}"

def run_mafft(sequences, output_format="clustal"):
    """
    Run MAFFT on a set of sequences for multiple sequence alignment.
    
    Args:
        sequences (dict): Dictionary with sequence names as keys and sequences as values
        output_format (str): Output format ('clustal' for clustal format, 'default' for fasta format)
        
    Returns:
        str: The alignment result
    """
    # Create a temporary FASTA file with input sequences
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as temp_in:
        for name, seq in sequences.items():
            temp_in.write(f">{name}\n{seq}\n")
        temp_in_name = temp_in.name
    
    # Create a temporary output file
    temp_out_handle, temp_out_name = tempfile.mkstemp(suffix='.aln')
    os.close(temp_out_handle)
    
    try:
        # Configure MAFFT command
        cmd = ["mafft", "--quiet"]
        
        # Add output format if specified
        if output_format == "clustal":
            cmd.append("--clustalout")
        
        cmd.extend([temp_in_name])
        
        # Log the MAFFT command being executed
        print(f"Running MAFFT command: {' '.join(cmd)}")
        
        # Execute MAFFT and capture output
        process = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            check=True
        )
        
        # Store the alignment result from stdout
        alignment_result = process.stdout
            
    except subprocess.CalledProcessError as e:
        error_message = f"MAFFT execution failed with error: {e.stderr}"
        print(error_message)
        return error_message
    except FileNotFoundError:
        error_message = "MAFFT command not found. Please ensure MAFFT is installed and in your system's PATH."
        print(error_message)
        return error_message
    finally:
        # Clean up temporary files
        try:
            os.remove(temp_in_name)
            os.remove(temp_out_name) 
        except Exception as e:
            print(f"Warning: Could not remove temporary files: {e}")
        
    return alignment_result.strip()

class InteractiveInterpreter(code.InteractiveInterpreter):
    """
    Custom interactive interpreter for executing Python code.
    Captures stdout and stderr for return to the UI.
    """
    def __init__(self):
        super().__init__()
        self.stdout = io.StringIO()
        self.stderr = io.StringIO()

    def execute(self, code_str):
        """
        Execute a string of Python code and capture the output.
        
        Args:
            code_str (str): The Python code to execute.
            
        Returns:
            str: The result of execution or error message.
        """
        # Redirect stdout/stderr to our string buffers
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        
        try:
            # Execute the code
            self.runcode(code_str)
            
            # Capture output
            output = self.stdout.getvalue()
            error = self.stderr.getvalue()
            
            # Return appropriate message
            if error:
                return f"Error: {error.strip()}"
            return output.strip() if output else "No result returned"
        finally:
            # Always restore stdout/stderr and reset buffers
            sys.stdout = original_stdout
            sys.stderr = original_stderr
            self.stdout = io.StringIO()
            self.stderr = io.StringIO()

# Initialize the Python code interpreter with access to global functions and modules
interpreter_globals = globals().copy()
interpreter_globals['__name__'] = '__main__'  # Allow code to run as if in main module

# Create the interpreter and update its namespace
interpreter = InteractiveInterpreter()
interpreter.locals.update(interpreter_globals)

def execute_python_code(code_str):
    """
    Execute Python code using the InteractiveInterpreter.
    
    Args:
        code_str (str): The Python code to execute.
        
    Returns:
        str: The output of the execution.
    """
    if not code_str.strip():
        return "No code provided to execute."
    
    print(f"Executing code:\n{code_str[:200]}{'...' if len(code_str) > 200 else ''}")
    return interpreter.execute(code_str)

def change_markdown_image(text: str):
    """
    Modify markdown image syntax for Gradio display.
    
    Args:
        text (str): Original markdown text
        
    Returns:
        str: Modified markdown text with adjusted image paths
    """
    if not text or not isinstance(text, str):
        return text
    
    # Convert single-quoted image paths to Gradio's /file format
    modified_text = re.sub(r"!\[(.*?)\]\(\'(.*?)\'\)", r"![\1](/file=\2)", text)
    return modified_text

def gradio_launch():
    """
    Launch the Gradio interface for the EPIC chatbot.
    Sets up the UI components and event handlers.
    
    Returns:
        gr.Blocks: The configured Gradio web interface
    """
    with gr.Blocks(theme=gr.themes.Soft()) as demo:
        # Setup UI components
        gr.Markdown("""
                    # EPIC - Enhanced Platform for Innovative Chemistry
                    ## A specialized assistant for quantum dot barcode (QDB) assay development
                    
                    Design and optimize probes for nucleic acid detection with expert guidance
                    """)
        chatbot = gr.Chatbot(height=820, label="Conversation")
        msg = gr.Textbox(label='User Input', placeholder="Ask about QDB assay design or type !code followed by Python code to execute...")
        clear = gr.Button("Clear Chat")

        conversation_history = [{"role": "system", "content": """You are a helpful assistant with expert knowledge in nucleic acid sequence design. 
                                At the beginning of each interaction, It will then request from the user the design objective and target if applicable. 
                            You will then ask users to specify their specific design criteria for their nucleic acid designs, including preferences for 
                            sequence length, GC content, ending nucleotides, secondary structures, and avoidance of repetitive sequences, emphasizing that 
                            no action will be attempted without these critical inputs from the user. 
                            Once the design criteria and the design objective have been provided, you will outline the planned actions based on these inputs 
                            before proceeding with the design process. 
                                 You will generate a Python script to execute the planned actions. 
                                 You should use a relevant library such as biopython to calculate relevant parameters such as melting temperature.  
                                 You should always use generate a Python script and use external software libraries such as biopython and primer3-py for producing deterministic results. 
                                 Use gc_fraction from Biopython to calculate GC content.
                                 Use Tm_NN from Biopython to calculate melting temperature for QDB assays. 
                                 Note that when using gc_fraction, you need to multiply by 100 to get percentage values. The correct sequence access pattern is alignment[i].seq[start:end] NOT alignment[i, j].
                                 Use the "Bio.Seq module" to handle nucleic acid sequence related operations. 
                                 Use primer3-py to design PCR primers.   
                                 The reverse complement of a sequence is the sequence that is complementary to the original sequence, but with the bases in reverse order. It is generally shown in 5' to 3' direction.
                                                                
                                You should follow the default design criteria and guidelines for designing nucleic acid quantum dot barcode (QDB) assays unless otherwise specified by the user:
                                 1. Count and recheck the number of nucleotides in the target sequence. 
                                 2. Split the target sequence into 2 equal-length halves. Check that both halves have an equal number of nucleotides. Then double-check the second half has the same length as the first half before you proceed. If the target sequence length is odd-numbered, leave one nucleotide in the middle as a spacer. 
                                 3. Reverse complement each half of the split target sequence and double check that they are correct
                                 4. To generate the reverse complement while considering directionality, remember that the 5' end of the new sequence will correspond to the 3' end of the original sequence, and vice versa. So if the original sequence is 5'-ATGC-3', the reverse complement will be 5'-GCAT-3'. If the original sequence is 3'-ATGC-5', the reverse complement will be 5'-TACG-3'. 
                                 5. Choose the second half of the reverse complement to be the capture probe and the first half to be the reporter probe
                                 6. Check the GC content of the capture and reporter probes, they should be between 40% and 60%, ideally 50% .
                                 7. Check the melting temperature of the capture and reporter probes. They should be between 55°C and 72°C, ideally 65°C.
                                 8. Lastly add an amine terminus denoted by "/5AmMC6/" to the 5' end of the capture probe, then add the Cy5 fluorophore denoted by "/3Cy5Sp/" to the 3' end of the reporter probe. 

                                 You should follow the default design criteria and guidelines for designing PCR primers unless otherwise specified by the user:
                                 - In general, a length of 18–30 nucleotides for primers is good.
                                 - Try to make the melting temperature (Tm) of the primers between 50°C and 65°C, and within 5°C of each other.
                                 - If the Tm of your primer is very low, try to find a sequence with more GC content, or extend the length of the primer a little.
                                 - Aim for the GC content to be between 40% and 60%, with the 3' of a primer ending in C or G to promote binding.
                                 - Typically, 3 to 4 nucleotides are added 5’ of the restriction enzyme site in the primer to allow for efficient cutting.
                                 - Try to avoid regions of secondary structure, and have a balanced distribution of GC-rich and AT-rich domains.
                                 - Try to avoid runs of 4 or more of one base, or dinucleotide repeats (for example, ACCCC or ATATATAT).
                                 - Avoid intra-primer homology (more than 3 bases that complement within the primer) or inter-primer homology (forward and reverse primers having complementary sequences). These circumstances can lead to self-dimers or primer-dimers instead of annealing to the desired DNA sequences.
                                 - The amplicon should start with the forward primer sequence, the start of the amplicon's reverse complement should be the reverse primer. 
                                 - Use the "design_primers" function to design PCR primers with Primer3-py.
                                 - Here is the example usage of primer3's 'design_primers' function usage: bindings.design_primers(    seq_args={        'SEQUENCE_ID': 'MH1000',        'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTT'                             'AGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCA'                             'ACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACG'                             'CACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAG'                             'TTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAA'                             'TGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAA'                             'ATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAA'                             'TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA'                             'GATGTTTCCCTCTAGTAG',        'SEQUENCE_INCLUDED_REGION': [36,342]    },    global_args={        'PRIMER_OPT_SIZE': 20,        'PRIMER_PICK_INTERNAL_OLIGO': 1,        'PRIMER_INTERNAL_MAX_SELF_END': 8,        'PRIMER_MIN_SIZE': 18,        'PRIMER_MAX_SIZE': 35,        'PRIMER_OPT_TM': 60.0,        'PRIMER_MIN_TM': 55.0,        'PRIMER_MAX_TM': 72.0,        'PRIMER_MIN_GC': 20.0,        'PRIMER_MAX_GC': 80.0,        'PRIMER_MAX_POLY_X': 100,        'PRIMER_INTERNAL_MAX_POLY_X': 100,        'PRIMER_SALT_MONOVALENT': 50.0,        'PRIMER_DNA_CONC': 50.0,        'PRIMER_MAX_NS_ACCEPTED': 0,        'PRIMER_MAX_SELF_ANY': 12,        'PRIMER_MAX_SELF_END': 8,        'PRIMER_PAIR_MAX_COMPL_ANY': 12,        'PRIMER_PAIR_MAX_COMPL_END': 8,        'PRIMER_PRODUCT_SIZE_RANGE': [            [75,100],[100,125],[125,150],            [150,175],[175,200],[200,225]        ],    })

                                
                                
                                You should follow the following guidelines when generating Python scripts with file operations:
                                    1. Do not use a main() function and the 'if __name__ == "__main__":' block.
                                    2. Use absolute paths: os.path.join(os.path.dirname(os.path.abspath(__file__)), "filename")
                                    3. Check directory existence: os.makedirs(directory, exist_ok=True)
                              
                                Only use MAFFT for multiple sequence alignment. 

                                MAFFT INTEGRATION:
                                 CRITICAL: A function called 'run_mafft' is already defined in this environment. 
                                 
                                 REQUIRED IMPORT FOR MAFFT:
                                 from docker_chat import run_mafft
                                 
                                 The run_mafft function takes two parameters:
                                    - sequences: A dictionary where keys are sequence names and values are the sequence strings
                                    - output_format: Either "clustal" for CLUSTAL format or "default" for FASTA format (default is "clustal")
                                 
                                 
                                 """}]

        def bot(history):
            """
            Handle messages from the user and generate bot responses.
            Streams LLM output and automatically executes any Python code.
            
            Args:
                history (list): Chat history
                
            Yields:
                list: Updated chat history for Gradio to render incrementally
            """
            # Get the user's message from history
            user_message = history[-1][0]
            user_message_str = str(user_message) if user_message is not None else ""
            
            # Add user message to conversation history for the LLM
            conversation_history.append({"role": "user", "content": user_message_str})

            # Handle direct code execution (!code prefix)
            if user_message_str.startswith("!code"):
                code_str = user_message_str[len("!code"):].strip()
                result = execute_python_code(code_str)
                history[-1][1] = f"Executing direct command:\n```python\n{code_str}\n```\n\n**Result:**\n```\n{result}\n```"
                yield history 
                return  # End generator here for !code command
                
            # Handle normal LLM interaction
            else:
                # Initialize bot's response for streaming
                history[-1][1] = "" 
                full_response = ""
                
                # Stream the response from GPT-4
                for chunk in chat_with_gpt4(conversation_history):
                    if chunk: 
                        history[-1][1] += chunk
                        full_response += chunk
                        yield history  # Yield updated history for each chunk

                # Save the assistant's response to conversation history if we got something back
                if full_response:
                    conversation_history.append({"role": "assistant", "content": full_response})

                # Auto-execute any Python code blocks in the response
                code_blocks = re.findall(r"```python\n(.*?)\n```", full_response, re.DOTALL)
                if code_blocks:
                    code_str_from_llm = code_blocks[0].strip()  # Take the first code block
                    if code_str_from_llm:  # Ensure extracted code is not empty
                        # Show executing message
                        history.append([None, "Executing generated code..."])
                        yield history 

                        # Execute the code and show results
                        result = execute_python_code(code_str_from_llm)
                        history[-1][1] = f"Execution Result:\n```\n{result}\n```"
                        yield history

        def user(user_message, history):
            """
            Process user input and update chat history.
            
            Args:
                user_message (str): User's input message
                history (list): Current chat history
                
            Returns:
                tuple: Empty string (to clear input box) and updated history
            """
            if not user_message:
                return "", history  # Don't add empty messages to history
            
            return "", history + [[user_message, None]]

        # Setup event handlers
        msg.submit(user, [msg, chatbot], [msg, chatbot], queue=False).then(
            bot, chatbot, chatbot
        )
        clear.click(lambda: [], None, chatbot, queue=False)  # Return empty list to clear chat

    # Configure and launch the Gradio interface
    demo.queue()
    demo.launch(
        share=True,          # For temporary sharing
        server_name="0.0.0.0",  # Binds to all network interfaces
        server_port=7860    # Standard Gradio port
    )
    return demo

if __name__ == "__main__":
    """
    Main entry point for the application.
    Checks prerequisites and launches the Gradio interface.
    """
    # Check for OpenAI API key
    if not os.getenv("OPENAI_API_KEY"):
        print("Error: OpenAI API key not found. Please set the OPENAI_API_KEY environment variable.")
        print("You can create a .env file with OPENAI_API_KEY=your_api_key")
        sys.exit(1)
    
    # Check for data directory
    data_dir = Path(__file__).parent / "data"
    if not data_dir.exists():
        print(f"Creating data directory at {data_dir}")
        data_dir.mkdir(exist_ok=True)
    
    # Check for MAFFT installation
    try:
        subprocess.run(["mafft", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("MAFFT installation detected.")
    except FileNotFoundError:
        print("Warning: MAFFT not found. Please install MAFFT for sequence alignment functionality.")
        print("Installation instructions: https://mafft.cbrc.jp/alignment/software/")
    
    print("Starting EPIC chatbot...")
    print("Access the web interface at http://localhost:7860")
    gradio_launch()
