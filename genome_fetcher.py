import os
from pathlib import Path
from typing import List, Union
from Bio import Entrez, SeqIO
import requests
from Bio.SeqRecord import SeqRecord
from .utils.logger import setup_logger
from .utils.exceptions import GenomeFetchError

logger = setup_logger(__name__)

class GenomeFetcher:
    """Fetch genomes from NCBI or load local files."""
    
    def __init__(self, email: str = "your@email.com"):
        self.email = email
        Entrez.email = email

    def fetch_ncbi(
        self, 
        search_term: str, 
        db: str = "nucleotide", 
        limit: int = 10
    ) -> List[SeqRecord]:
        """Fetch genomes from NCBI using a search term."""
        try:
            handle = Entrez.esearch(db=db, term=search_term, retmax=limit)
            id_list = Entrez.read(handle)["IdList"]
            handle = Entrez.efetch(db=db, id=id_list, rettype="gb", retmode="text")
            return list(SeqIO.parse(handle, "genbank"))
        except Exception as e:
            logger.error(f"NCBI fetch failed: {e}")
            raise GenomeFetchError(f"NCBI API error: {e}")

    def load_local(
        self, 
        dir_path: Union[str, Path], 
        file_format: str = "fasta"
    ) -> List[SeqRecord]:
        """Load genomes from a local directory."""
        try:
            dir_path = Path(dir_path)
            if not dir_path.exists():
                raise FileNotFoundError(f"Directory {dir_path} not found")
            return [
                SeqIO.read(f, file_format) 
                for f in dir_path.glob(f"*.{file_format}")
            ]
        except Exception as e:
            logger.error(f"Local load failed: {e}")
            raise GenomeFetchError(f"Local load error: {e}")
