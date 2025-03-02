# test_api.py
from genome_fetcher import GenomeFetcher

fetcher = GenomeFetcher(
    email="test@example.com",
    api_key="test_key"
)
print("✅ Successfully initialized GenomeFetcher")
