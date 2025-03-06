# GuideX
crispr/cas13 software that identifies conserved regions, designs cas13 gRNAs, performs off-target analysis, and completes gRNA optimization.

## Requirements
- MUSCLE v5+ (`brew install muscle`)
- Python 3.10+
- Dependencies: `pip install -r requirements.txt`

## Usage
```bash
export NCBI_API_KEY_2025="your_ncbi_key"
python3 test_guidex.py
