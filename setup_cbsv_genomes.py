import os
import subprocess
from pathlib import Path

# Configuration
ACCESSIONS = [
    "NC_012698.2", "GU563327.1", "FN434436.1", "MK955888.1",
    "OQ335847.1", "KY563367.1", "MG019915.1", "MK103393.1",
    "MW961170.2", "MZ362877.1", "GU563326.1", "KY290025.1",
    "LT577538.1", "FN434437.1"
]

GENOME_DIR = Path("genomes/CBSV")
INDEX_DIR = GENOME_DIR / "bowtie_indexes"
INDEX_DIR.mkdir(parents=True, exist_ok=True)

def install_tools():
    """Install required NCBI utilities"""
    try:
        subprocess.run(["efetch", "-h"], check=True,
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except:
        print("⏳ Installing Entrez Direct...")
        subprocess.run(
            'sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"',
            shell=True, check=True
        )
    
    # Add edirect to PATH
    subprocess.run('export PATH="${PATH}:${HOME}/edirect"', shell=True)

def download_genome(accession: str, idx: int):
    """Download sequence using efetch"""
    output_file = GENOME_DIR / f"CBSV_{idx+1}.fna"
    
    if output_file.exists():
        print(f"✓ Genome {accession} already downloaded")
        return output_file

    print(f"⏬ Downloading {accession}...")
    try:
        with open(output_file, "w") as f:
            subprocess.run([
                "efetch", 
                "-db", "nuccore",
                "-id", accession,
                "-format", "fasta"
            ], check=True, stdout=f, text=True)
            
        return output_file
    except Exception as e:
        print(f"❌ Failed to download {accession}: {str(e)}")
        return None

def build_index(fasta_path: Path, idx: int):
    """Build Bowtie2 index (same as before)"""
    # ... keep existing build_index function ...

if __name__ == "__main__":
    install_tools()
    GENOME_DIR.mkdir(parents=True, exist_ok=True)
    
    # Update PATH for current session
    os.environ["PATH"] += os.pathsep + os.path.expanduser("~/edirect")
    
    # Download and index genomes
    for idx, accession in enumerate(ACCESSIONS):
        genome_path = download_genome(accession, idx)
        if genome_path:
            build_index(genome_path, idx)

    print("\n🎉 Setup complete! Verify with:")
    print(f"ls -lh {GENOME_DIR}/*.fna")
    print(f"ls -lh {INDEX_DIR}/*.bt2")
