import os
from pathlib import Path
from typing import List, Union
import requests
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from .utils.logger import setup_logger
from .utils.exceptions import GenomeFetchError
import time
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = setup_logger(__name__)

class GenomeFetcher:
    """Fetch genomes from NCBI Datasets API v2 with modern features"""
    
    API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
    RATE_LIMIT = 5  # requests per second
    RETRY_STRATEGY = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504]
    )
    
    def __init__(self, email: str, api_key: str = None):
        self.email = email
        self.api_key = api_key
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=self.RETRY_STRATEGY))

    def fetch_ncbi(
        self, 
        search_term: str,
        limit: int = 5,
        exclude_atypical: bool = True
    ) -> List[SeqRecord]:
        """Fetch genomes using NCBI Datasets API v2"""
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
                self.RATE_LIMIT = 10  # Increased limit with API key

            # Get dehydrated dataset
            resp = self._rate_limited_request(
                f"{self.API_BASE}/genome/search",
                params=params
            )
            data = resp.json()
            
            # Rehydrate sequences
            return self._process_dataset(data["data"]["reports"])
            
        except Exception as e:
            logger.error(f"API v2 fetch failed: {e}")
            raise GenomeFetchError(f"NCBI API v2 error: {e}")

    def _rate_limited_request(self, url, **kwargs):
        """Handle NCBI's rate limits with exponential backoff"""
        delay = 1 / self.RATE_LIMIT
        time.sleep(delay)
        response = self.session.get(url, timeout=15, **kwargs)
        response.raise_for_status()
        return response

    def _process_dataset(self, reports: list) -> List[SeqRecord]:
        """Convert API response to SeqRecords with modern filtering"""
        genomes = []
        for report in reports:
            try:
                if self._is_valid_genome(report):
                    seq = Seq(report["sequence"]["sequence"])
                    genomes.append(SeqRecord(
                        seq,
                        id=report["accession"],
                        description=f"{report['organism']['name']} | {report['length']}bp"
                    ))
            except KeyError as e:
                logger.warning(f"Skipping invalid genome report: {e}")
        return genomes

    def _is_valid_genome(self, report: dict) -> bool:
        """Apply modern validation criteria"""
        return (
            report["length"] >= 1000 and
            not report.get("is_atypical", False) and
            report["sequence"]["sequence"] is not None
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
