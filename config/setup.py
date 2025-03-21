from setuptools import setup, find_packages

setup(
    name="guidex",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "scikit-bio>=0.5.6", 
        "numpy>=1.21.0",
        "pandas>=1.3.5",
        "plotly>=5.10.0",
        "requests>=2.26.0",
        "scipy>=1.7.0",
        "pyyaml>=6.0",
        "importlib-resources>=5.4.0",  # Added for config access
    ],
    package_data={
        'guidex': ['../config/cas13_subtypes.yaml']  # Include config
    },
    include_package_data=True,
    python_requires=">=3.8",   # Specify Python version
    entry_points={             # Optional but useful
        "console_scripts": [
            "guidex=guidex.cli:main",
        ],
    },
    author="Darsh Shrivastava",
    author_email="darsh.shri123@email.com",
    description="CRISPR guide RNA design system with conservation analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/guidex",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
