from typing import List, Optional
import requests
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
import time
import os

class GenomeFetcher:
    """Modern NCBI Datasets API v2 integration with full API key support"""
    
    API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
    RATE_LIMIT = 5  # Default requests/sec
    RETRY_STRATEGY = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504]
    )

    def __init__(self, email: str, api_key: Optional[str] = None):
        """
        Initialize with API credentials
        :param email: Required for NCBI compliance
        :param api_key: Optional for higher rate limits
        """
        self.email = email
        self.api_key = api_key
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=self.RETRY_STRATEGY))
        
        # Enhanced rate limiting with API key
        self.RATE_LIMIT = 10 if api_key else 5

    def fetch_ncbi(
        self, 
        search_term: str,
        limit: int = 5,
        exclude_atypical: bool = True
    ) -> List[SeqRecord]:
        """
        Fetch genomes using NCBI Datasets API v2
        :param search_term: e.g. "Influenza A virus[Organism]"
        :param limit: Maximum results to return
        :param exclude_atypical: Filter out problematic genomes
        """
        try:
            params = {
                "term": search_term,
                "limit": limit,
                "data_type": "genome",
                "include": "sequence",
                "exclude_atypical": str(exclude_atypical).lower()
            }
            
            if self.api_key:
                params["api_key"] = self.api_key

            response = self._rate_limited_request(
                f"{self.API_BASE}/genome/search",
                params=params
            )
            return self._process_response(response.json())
            
        except Exception as e:
            raise RuntimeError(f"API v2 fetch failed: {str(e)}")

    def _rate_limited_request(self, url, params=None):
        """Handle NCBI rate limits with exponential backoff"""
        time.sleep(1/self.RATE_LIMIT)  # Enforce rate limit
        response = self.session.get(
            url,
            params=params,
            headers={"User-Agent": f"GuideX/{os.getenv('VERSION', '1.0')}"},
            timeout=15
        )
        response.raise_for_status()
        return response

    def _process_response(self, data: dict) -> List[SeqRecord]:
        """Convert API response to validated SeqRecords"""
        genomes = []
        for item in data.get("data", {}).get("reports", []):
            try:
                if self._is_valid_genome(item):
                    seq = self._create_sequence(item)
                    genomes.append(seq)
            except KeyError as e:
                continue
        return genomes

    def _is_valid_genome(self, genome_data: dict) -> bool:
        """Apply NCBI's quality criteria"""
        return (
            genome_data.get("length", 0) >= 1000 and
            not genome_data.get("is_atypical", False) and
            "sequence" in genome_data and
            "sequence" in genome_data["sequence"]
        )

    def _create_sequence(self, genome_data: dict) -> SeqRecord:
        """Create standardized SeqRecord from API data"""
        return SeqRecord(
            Seq(genome_data["sequence"]["sequence"]),
            id=genome_data.get("accession", "UNKNOWN"),
            description=f"{genome_data.get('organism', {}).get('name', '')} | "
                       f"{genome_data.get('length', 0)}bp"
        )

def load_dehydrated(self, zip_path: Path) -> List[SeqRecord]:
    """Load and rehydrate genomes from NCBI dehydrated package"""
    from zipfile import ZipFile
    import tempfile
    import json
    
    genomes = []
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Extract package
            with ZipFile(zip_path, 'r') as z:
                z.extractall(tmpdir)
                
            # Read dataset catalog
            catalog_path = Path(tmpdir) / "ncbi_dataset" / "dataset_catalog.json"
            with open(catalog_path) as f:
                catalog = json.load(f)
                
            # Process sequences
            for assembly in catalog["assemblies"]:
                seq_file = Path(tmpdir) / assembly["sequenceFile"]
                if seq_file.exists():
                    with open(seq_file) as sf:
                        for record in SeqIO.parse(sf, "fasta"):
                            genomes.append(record)
                            if len(genomes) >= 100:  # Safety limit
                                return genomes
                                
            return genomes
            
    except Exception as e:
        logger.error(f"Dehydrated load failed: {e}")
        raise GenomeFetchError(f"Dehydrated package error: {e}")
