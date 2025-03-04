import requests
from typing import List, Optional, Dict, Union
from pathlib import Path
import subprocess
import json
import tempfile
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import re
import os
import logging
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
import time

class GenomeFetcher:
    """Modern NCBI Genome Fetcher using CLI v16+ and API v2 (2025 standards)"""
    
    CLI_PATH = "/usr/local/ncbi/datasets"  # <-- PROPER CLASS LEVEL
    ACCESSION_REGEX = r"^GC[AFN]_[0-9]{11}\.\d$"
    RATE_LIMIT = 8
    RETRY_STRATEGY = Retry(
        total=5,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504]
    )

    def __init__(self, api_key: Optional[str] = None):  # Line 18
        # PROPERLY INDENTED INIT CODE
        self.api_key = api_key or os.getenv("NCBI_API_KEY_2025")
        self._verify_cli()
        self.logger = logging.getLogger(__name__)
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=self.RETRY_STRATEGY))
        self._configure_cli()
            
    def _verify_cli(self):
        """Validate CLI v16+ installation"""
        try:
            result = subprocess.run(
                ["datasets", "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            if "16." not in result.stdout:
                raise RuntimeError("Requires datasets CLI v16+")
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError(
                "NCBI datasets CLI not found. Install with: "
                "brew install ncbi-datasets"
            )
    def _configure_cli(self):
        """Configure CLI authentication"""
        if self.api_key:
            try:
                subprocess.run(
                    [self.CLI_PATH, "config", "set", "api-key", self.api_key],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )
            except subprocess.CalledProcessError as e:
                self.logger.warning("Failed to configure CLI authentication")

    def fetch_genomes(
        self,
        target: Union[str, List[str]],
        genome_type: str = "reference"
    ) -> List[SeqRecord]:
        """Unified fetch method for 2025 standards"""
        if isinstance(target, list) and all(re.match(self.ACCESSION_REGEX, t) for t in target):
            return self._fetch_by_accessions(target)
        elif isinstance(target, str):
            return self._fetch_by_taxonomy(target, genome_type)
        raise ValueError("Invalid target type - must be accession list or organism name")

    def _fetch_by_taxonomy(self, organism: str, genome_type: str) -> List[SeqRecord]:
        """Fetch genomes using 2025 taxonomic search"""
        cmd = [
            self.CLI_PATH, "summary", "genome", "taxon",
            organism,
            "--assembly-level", "chromosome",
            "--assembly-source", genome_type,
            "--annotated", 
            "--limit", "10",
            "--as-json-lines"
        ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=45
            )
            return self._process_jsonl(result.stdout)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Taxonomy fetch failed: {e.stderr}")
            return self._fallback_strategy()

    def _fetch_by_accessions(self, accessions: List[str]) -> List[SeqRecord]:
        """Batch fetch using 2025 accession format"""
        with tempfile.TemporaryDirectory() as tmpdir:
            accession_file = Path(tmpdir) / "accessions.txt"
            accession_file.write_text("\n".join(accessions))
            
            cmd = [
                self.CLI_PATH, "download", "genome",
                "accession", "--inputfile", str(accession_file),
                "--include", "genome,cds,rna",
                "--filename", f"{tmpdir}/dataset.zip"
            ]
            
            try:
                subprocess.run(cmd, check=True, timeout=120)
                return self._extract_genomes(f"{tmpdir}/dataset.zip")
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"Accession download failed: {e.stderr}")
                
    def _rate_limited_request(self, url, params=None, headers=None):  # Add headers parameter
        """Handle NCBI rate limits with exponential backoff"""
        time.sleep(1/self.RATE_LIMIT)
        response = self.session.get(
            url,
            params=params,
            headers=headers,  # Pass headers here
            timeout=15
        )
        response.raise_for_status()
        return response
    
    def _process_jsonl(self, jsonl_data: str) -> List[SeqRecord]:
        """Process 2025 JSON Lines format with quality control"""
        genomes = []
        for line in jsonl_data.splitlines():
            data = json.loads(line)
            if self._meets_quality_standards(data):
                genomes.append(self._create_record(data))
        return genomes

    def _meets_quality_standards(self, data: Dict) -> bool:
        """2025 genome quality criteria"""
        return (
            data.get("assembly", {}).get("assembly_level") == "chromosome" and
            data.get("annotation", {}).get("quality") == "gold" and
            data.get("assembly", {}).get("contig_n50") >= 100000
        )

    def _create_record(self, data: Dict) -> SeqRecord:
        """Create annotated SeqRecord from 2025 schema"""
        features = [
            f"{feat['type']}:{feat['location']}" 
            for feat in data.get("features", [])
        ]
        
        return SeqRecord(
            Seq(data["assembly"]["sequence"]),
            id=data["assembly"]["accession"],
            description=f"{data['organism']['name']} | {data['assembly']['name']} | " +
                       f"Features: {', '.join(features)}"
        )

    def _extract_genomes(self, zip_path: str) -> List[SeqRecord]:
        """Extract and validate genomes from 2025 package format"""
        genomes = []
        with tempfile.TemporaryDirectory() as tmpdir:
            subprocess.run(
                ["unzip", zip_path, "-d", tmpdir],
                check=True,
                stdout=subprocess.DEVNULL
            )
            
            seq_dir = Path(tmpdir) / "ncbi_dataset" / "data"
            for seq_file in seq_dir.glob("**/*.fna"):
                for record in SeqIO.parse(seq_file, "fasta"):
                    if len(record.seq) >= 1000:
                        genomes.append(record)
        return genomes

    def _fallback_strategy(self):
        """2025 resilient fallback to cached genomes"""
        cache_path = Path("~/.guidex/genome_cache").expanduser()
        if cache_path.exists():
            self.logger.info("Using cached genome backup")
            return list(SeqIO.parse(cache_path / "default_genomes.fna", "fasta"))
        raise RuntimeError("No genomes available and fallback cache missing")

    def preload_cache(self, genomes: List[SeqRecord], organism: str):
        """Cache genomes with organism-specific organization"""
        cache_dir = Path("~/.guidex/genome_cache").expanduser()
        organism_dir = cache_dir / organism.replace(" ", "_")
        organism_dir.mkdir(parents=True, exist_ok=True)
    
        cache_file = organism_dir / "genomes.fna"
    
        with open(cache_file, "w") as f:
            SeqIO.write(genomes, f, "fasta-2line")
    
        self.logger.info(f"Cached {len(genomes)} genomes for {organism} at {cache_file}")

    def load_cached_genomes(self, organism: str) -> List[SeqRecord]:
        """Load organism-specific genomes from cache"""
        cache_dir = Path("~/.guidex/genome_cache").expanduser()
        organism_dir = cache_dir / organism.replace(" ", "_")
        cache_file = organism_dir / "genomes.fna"
    
        if cache_file.exists():
            self.logger.debug(f"Loading cached genomes from {cache_file}")
            return list(SeqIO.parse(cache_file, "fasta"))
        raise FileNotFoundError(f"No cached genomes for {organism}")

    @property
    def status(self) -> Dict:
        """Get API/cli status"""
        result = subprocess.run(
            [self.CLI_PATH, "version"],
            capture_output=True,
            text=True
        )
        return {
            "cli_version": result.stdout.strip(),
            "api_status": "active" if self.api_key else "unauthenticated",
            "rate_limit": self.RATE_LIMIT
        }

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
