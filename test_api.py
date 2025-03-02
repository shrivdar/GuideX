# In test_api.py
import os
from genome_fetcher import GenomeFetcher

# Get key from environment
API_KEY = os.getenv("NCBI_API_KEY") 

fetcher = GenomeFetcher(
    email="darsh.shri123@gmail.com",
    api_key=API_KEY  # Pass from environment
)
