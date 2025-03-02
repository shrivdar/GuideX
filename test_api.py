import os
from genome_fetcher import GenomeFetcher

def test_constructor():
    # Test without API key
    fetcher = GenomeFetcher(email="test@example.com")
    print("✅ Constructor works without API key")
    
    # Test with API key
    fetcher = GenomeFetcher(
        email="test@example.com",
        api_key="dummy_key_123"  # Test value
    )
    print("✅ Constructor works with API key")

if __name__ == "__main__":
    test_constructor()
