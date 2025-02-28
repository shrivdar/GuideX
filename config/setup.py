from setuptools import setup, find_packages

setup(
    name="guidex",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "scikit-bio>=0.5.6",
        "pandas>=1.3.5",
        "plotly>=5.10.0",
        "pyyaml>=6.0",
    ],
    entry_points={
        "console_scripts": [
            "guidex=guidex.cli:main",  # Add a CLI entry point later
        ],
    },
)
