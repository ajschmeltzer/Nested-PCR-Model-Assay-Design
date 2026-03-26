PCR Assay Designer

Automated nested PCR assay design pipeline built with Python and Primer3.
The pipeline downloads genomes (or specific genome regions) from NCBI and designs two-level nested PCR assays with annotated outputs.
-----------------------------------------------------------------------
Features:
	-Automatic genome download from NCBI
	-Support for user-defined genome regions (e.g., NC_000913.3:3734-5020)
	-Local genome/region caching to avoid repeated downloads
	-Random genome region sampling (for full genomes)
	-Automated outer and inner nested PCR primer design
	-Configurable Primer3 parameters via config.yaml
	-Multiple assays per organism (full genome mode)
	-Annotated ApE plasmid files for visualization
	-CSV and TXT primer reports
-----------------------------------------------------------------------
Installation

Clone the repository:
git clone https://github.com/ajschmeltzer/Nested-PCR-Model-Assay-Design
cd Nested-PCR-Model-Assay-Design

Install dependencies:
pip install -r requirements.txt

Dependencies include:
	Biopython
	primer3-py
	PyYAML
-----------------------------------------------------------------------
Usage

Run the pipeline from the repository root:
python run_pipeline.py

The pipeline will:
	-Download genomes or genome regions if not already present
	-If full genome: sample random genome regions
	-If region specified: use that exact region
	-Design nested PCR assays
	-Export results
-----------------------------------------------------------------------
Notes:

Regions are fetched directly from NCBI using seq_start and seq_stop
When a region is specified:
	-Only one assay is generated
	-No random sampling is performed
	-The full region is used for primer design
-----------------------------------------------------------------------
Project Structure

Nested-PCR-Model-Assay-Design
│
├── config.yaml
├── run_pipeline.py
├── requirements.txt
├── README.txt
│
├── src
│   └── assay_pipeline.py
│
└── Designed Assays
	├── Genomes
	└── Assays by Organism
-----------------------------------------------------------------------
File and Folder Descriptions

config.yaml
Pipeline configuration including organisms and primer settings.

run_pipeline.py
Entry point for running the pipeline.

src/assay_pipeline.py
Core primer design pipeline.

Designed Assays/Genomes
Downloaded genomes and regions stored locally.

Designed Assays/Assays by Organism
Generated assay output folders.
-----------------------------------------------------------------------
Configuration

All pipeline settings are controlled through config.yaml.

Example config.yaml:

ORGANISMS:
	E_coli: NC_000913.3
	E_coli_thrC: NC_000913.3:3734-5020

pipeline:
	assays_per_org: 5	REMINDER: Only one assay will be designed for the specific E_coli_thrC: NC_000913.3:3734-5020 region
	region_length: 10000
	email: your_email@example.com

primer3:
	outer:
		PRIMER_OPT_SIZE: 20
		PRIMER_MIN_SIZE: 18
		PRIMER_MAX_SIZE: 25
		PRIMER_PRODUCT_SIZE_RANGE: [[230, 400]]
	inner:
		PRIMER_PRODUCT_SIZE_RANGE: [[120, 200]]
-----------------------------------------------------------------------
Organisms

The ORGANISMS section is a dictionary of organism names and NCBI genome inputs.

Supported formats:

Full genome/chromosome:
Organism_name : NCBI_accession

	Example:
	E_coli : NC_000913.3
	Human_chr10: NC_000010.11

Specific genome region:
Organism_name : NCBI_accession

	Example:
	E_coli_thrC : NC_000913.3:3734-5020
-----------------------------------------------------------------------
Output

Generated assays are stored in:
Designed Assays/Assays by Organism/

Each assay folder contains:
	info.txt
	*_primers.txt
	*_primers.csv
	*.ape
	*_feature_library.txt
-----------------------------------------------------------------------
File Descriptions

*_primers.txt
Human-readable primer report.

*_primers.csv
Structured primer data.

*.ape
Annotated ApE sequence file showing primer binding sites.

*_feature_library.txt
Feature library for importing primers into ApE.

info.txt
Genome location and metadata for the assay.
Indicates whether the assay was generated from:
- Random genome region
- User-specified region
-----------------------------------------------------------------------
Genome Storage

Downloaded genomes and regions are stored locally in:
Designed Assays/Genomes

Files are cached and reused for future runs.

Region files are saved as:
.fasta
-----------------------------------------------------------------------
Requirements

-Python 3.9 or newer
-Internet connection for initial genome download
-installation of packages included in the requirements.txt file
-----------------------------------------------------------------------
Future Improvements

-BLAT or BLAST specificity checking
-Masking of low-complexity or non-conserved genome regions
-Improved error handling for failed NCBI fetches
-----------------------------------------------------------------------
License

MIT License
