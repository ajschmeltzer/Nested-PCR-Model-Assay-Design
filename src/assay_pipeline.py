import os
import re
import random
from pathlib import Path
from Bio import Entrez, SeqIO
import primer3
import csv
import yaml


# -----------------------------
# Base folders
# -----------------------------
# Automatically set the base folder to the root of the repository
BASE_FOLDER = Path(__file__).parent.parent.resolve()
GENOME_FOLDER = BASE_FOLDER / "Designed Assays" / "Genomes"
OUTPUT_FOLDER = BASE_FOLDER / "Designed Assays" / "Assays by Organism"

# Make sure directories exist
GENOME_FOLDER.mkdir(parents=True, exist_ok=True)
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

#Set base folder
def get_base_folder():
    return Path(os.environ.get("PCR_ASSAY_ROOT", BASE_FOLDER))

# -----------------------------
# YAML Information
# -----------------------------
# Load config.yaml
# CONFIG_FILE = BASE_FOLDER / "config.yaml"
# with open(CONFIG_FILE) as cf:
#     config = yaml.safe_load(cf)

# # Pipeline settings
# ORGANISMS = config.get("organisms", {})
# PIPELINE_CONFIG = config.get("pipeline", {})
# ASSAYS_PER_ORG = PIPELINE_CONFIG.get("assays_per_org", 5)
# REGION_LENGTH = PIPELINE_CONFIG.get("region_length", 1000)
# EMAIL = PIPELINE_CONFIG.get("email", "your_email@example.com")

# # Primer3 settings
# PRIMER3_CONFIG = config.get("primer3", {})
# OUTER_SETTINGS = PRIMER3_CONFIG.get("outer", {})
# INNER_SETTINGS = PRIMER3_CONFIG.get("inner", {})

# -----------------------------
# Helper functions
# -----------------------------
def gc_content(seq):
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return 100 * gc_count / len(seq)

def revcomp(seq):
    complement = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(complement)[::-1]

# -----------------------------
# NCBI fetcher
# -----------------------------
class NCBISequenceFetcher:
    def __init__(self, email, genome_dir=None):
        Entrez.email = email
        if genome_dir is None:
            genome_dir = GENOME_FOLDER
        self.genome_dir = Path(genome_dir)
        self.genome_dir.mkdir(parents=True, exist_ok=True)

    def fetch_genome(self, accession):
        """Download full genome if missing, hard-masked (N's in repeats)."""
        genome_path = self.genome_dir / f"{accession}_hardmasked.fasta"
        if genome_path.exists():
            print(f"Genome already downloaded: {genome_path.name}")
            safe_id = re.sub(r'[:<>"/\\|?*]', '_', accession)
            return genome_path, safe_id

        print(f"Fetching full genome from NCBI (hard-masked): {accession} ...")
        handle = Entrez.efetch(
            db="nuccore",
            id=accession,
            rettype="fasta",
            retmode="text",
            mask="soft"  # still downloads lowercase/hard N's depending on NCBI
        )
        record = SeqIO.read(handle, "fasta")
        SeqIO.write(record, genome_path, "fasta")
        safe_id = re.sub(r'[:<>"/\\|?*]', '_', accession)
        print(f"Genome saved: {genome_path}")
        return genome_path, safe_id

# -----------------------------
# ApE file writer
# -----------------------------
def write_ape_file(assay_name, outer_amp_seq, primers, out_dir):

    primer_colors = {
        "outer_forward": "#FF0000",   # red
        "outer_reverse": "#FF0000",
        "inner_forward": "#00FFFF",   # cyan
        "inner_reverse": "#00FFFF"
    }

    ape_file = out_dir / f"{assay_name}.ape"
    seq_upper = outer_amp_seq.upper()

    with open(ape_file, "w") as f:

        f.write(f"LOCUS       {assay_name} {len(outer_amp_seq)} bp DNA linear\n")
        f.write(f"DEFINITION  Outer amplicon for {assay_name}\n")
        f.write("FEATURES             Location/Qualifiers\n")

        for name, primer_seq in primers.items():

            primer_seq = primer_seq.upper()
            rc = revcomp(primer_seq)
            color = primer_colors[name]

            pos = seq_upper.find(primer_seq)

            if pos != -1:

                start = pos + 1
                end = pos + len(primer_seq)

                f.write(f"     primer_bind     {start}..{end}\n")
                f.write(f"                     /label=\"{name}\"\n")
                f.write(f"                     /ApEinfo_fwdcolor={color}\n")
                f.write(f"                     /ApEinfo_revcolor={color}\n")

            else:

                pos = seq_upper.find(rc)

                if pos != -1:

                    start = pos + 1
                    end = pos + len(primer_seq)

                    f.write(f"     primer_bind     complement({start}..{end})\n")
                    f.write(f"                     /label=\"{name}\"\n")
                    f.write(f"                     /ApEinfo_fwdcolor={color}\n")
                    f.write(f"                     /ApEinfo_revcolor={color}\n")

        f.write("ORIGIN\n")

        for i in range(0, len(outer_amp_seq), 60):
            seq_chunk = outer_amp_seq[i:i+60].upper()
            f.write(f"{i+1:9} {seq_chunk}\n")

        f.write("//\n")

    return ape_file

# -----------------------------
# ApE feature library writer
# -----------------------------
def write_ape_feature_library(
    assay_name,
    outer_forward,
    outer_reverse,
    inner_forward,
    inner_reverse,
    assay_folder
):

    feature_file = assay_folder / f"{assay_name}_feature_library.txt"

    with open(feature_file, "w") as f:

        f.write(f"outer_forward\t{outer_forward}\tprimer_bind\tred\n")
        f.write(f"outer_reverse\t{outer_reverse}\tprimer_bind\tred\tred\n")
        f.write(f"inner_forward\t{inner_forward}\tprimer_bind\tcyan\n")
        f.write(f"inner_reverse\t{inner_reverse}\tprimer_bind\tcyan\tcyan\n")

    return feature_file

# -----------------------------
# Primer design function
# -----------------------------
def design_assay(sequence, safe_id, assay_folder, outer_settings, inner_settings):
    """
    Design nested PCR primers for a given sequence.
    outer_settings and inner_settings must be dictionaries from config.yaml.
    """

    # -----------------------------
    # Outer primers
    # -----------------------------
    outer_results = primer3.bindings.design_primers(
        {'SEQUENCE_ID': safe_id, 'SEQUENCE_TEMPLATE': sequence},
        outer_settings
    )

    outer_forward = outer_results['PRIMER_LEFT_0_SEQUENCE']
    outer_reverse = outer_results['PRIMER_RIGHT_0_SEQUENCE']
    outer_product_size = outer_results['PRIMER_PAIR_0_PRODUCT_SIZE']
    outer_left_start, outer_left_len = outer_results['PRIMER_LEFT_0']
    outer_right_start, outer_right_len = outer_results['PRIMER_RIGHT_0']
    outer_amplicon_seq = sequence[outer_left_start : outer_left_start + outer_product_size]

    outer_forward_tm = primer3.calc_tm(outer_forward)
    outer_reverse_tm = primer3.calc_tm(outer_reverse)
    outer_forward_gc = gc_content(outer_forward)
    outer_reverse_gc = gc_content(outer_reverse)
    outer_amplicon_tm = primer3.calc_tm(outer_amplicon_seq)

    # -----------------------------
    # Inner primers
    # -----------------------------
    inner_start = outer_left_start + outer_left_len + 5
    inner_end = outer_right_start - outer_right_len - 5
    inner_sequence = sequence[inner_start:inner_end]

    inner_results = primer3.bindings.design_primers(
        {'SEQUENCE_ID': safe_id + "_inner", 'SEQUENCE_TEMPLATE': inner_sequence},
        inner_settings
    )

    inner_forward = inner_results['PRIMER_LEFT_0_SEQUENCE']
    inner_reverse = inner_results['PRIMER_RIGHT_0_SEQUENCE']
    inner_product_size = inner_results['PRIMER_PAIR_0_PRODUCT_SIZE']
    inner_left_start, inner_left_len = inner_results['PRIMER_LEFT_0']
    inner_amplicon_seq = inner_sequence[inner_left_start : inner_left_start + inner_product_size]

    inner_forward_tm = primer3.calc_tm(inner_forward)
    inner_reverse_tm = primer3.calc_tm(inner_reverse)
    inner_forward_gc = gc_content(inner_forward)
    inner_reverse_gc = gc_content(inner_reverse)
    inner_amplicon_tm = primer3.calc_tm(inner_amplicon_seq)

    # ----------------------------- 
    # Save TXT
    # -----------------------------
    primer_file = assay_folder / f"{safe_id}_primers.txt"

    with open(primer_file, "w") as f:

        f.write("=== Outer Primers ===\n")
        f.write(f"Forward: {outer_forward.upper()}\n")
        f.write(f"Reverse: {outer_reverse.upper()}\n")
        f.write(f"Reverse (5'->3' revcomp): {revcomp(outer_reverse).upper()}\n")
        f.write(f"Product size: {outer_product_size}\n")
        f.write(f"Forward Tm: {outer_forward_tm:.2f} °C\n")
        f.write(f"Reverse Tm: {outer_reverse_tm:.2f} °C\n")
        f.write(f"Forward GC%: {outer_forward_gc:.2f}\n")
        f.write(f"Reverse GC%: {outer_reverse_gc:.2f}\n")
        f.write(f"Amplicon Tm (GC approx.): {outer_amplicon_tm:.2f} °C\n")
        f.write(f"Amplicon sequence: {outer_amplicon_seq.upper()}\n\n")

        f.write("=== Inner Primers ===\n")
        f.write(f"Forward: {inner_forward.upper()}\n")
        f.write(f"Reverse: {inner_reverse.upper()}\n")
        f.write(f"Reverse (5'->3' revcomp): {revcomp(inner_reverse).upper()}\n")
        f.write(f"Product size: {inner_product_size}\n")
        f.write(f"Forward Tm: {inner_forward_tm:.2f} °C\n")
        f.write(f"Reverse Tm: {inner_reverse_tm:.2f} °C\n")
        f.write(f"Forward GC%: {inner_forward_gc:.2f}\n")
        f.write(f"Reverse GC%: {inner_reverse_gc:.2f}\n")
        f.write(f"Amplicon Tm (GC approx.): {inner_amplicon_tm:.2f} °C\n")
        f.write(f"Amplicon sequence: {inner_amplicon_seq.upper()}\n")

    # -----------------------------
    # Save CSV
    # -----------------------------
    csv_file = assay_folder / f"{safe_id}_primers.csv"

    with open(csv_file, "w", newline="") as csvf:

        writer = csv.writer(csvf)

        writer.writerow([
            "Primer Type","Forward","Reverse","Reverse (revcomp)","Product Size",
            "Forward Tm","Reverse Tm","Forward GC","Reverse GC",
            "Amplicon Tm","Amplicon Sequence"
        ])

        writer.writerow([
            "Outer", outer_forward.upper(), outer_reverse.upper(), revcomp(outer_reverse).upper(),
            outer_product_size,
            f"{outer_forward_tm:.2f}", f"{outer_reverse_tm:.2f}",
            f"{outer_forward_gc:.2f}", f"{outer_reverse_gc:.2f}",
            f"{outer_amplicon_tm:.2f}", outer_amplicon_seq.upper()
        ])

        writer.writerow([
            "Inner", inner_forward.upper(), inner_reverse.upper(), revcomp(inner_reverse).upper(),
            inner_product_size,
            f"{inner_forward_tm:.2f}", f"{inner_reverse_tm:.2f}",
            f"{inner_forward_gc:.2f}", f"{inner_reverse_gc:.2f}",
            f"{inner_amplicon_tm:.2f}", inner_amplicon_seq.upper()
        ])

    # -----------------------------
    # Write ApE file (auto annotated)
    # -----------------------------
    ape_file = write_ape_file(
        safe_id,
        outer_amplicon_seq,
        {
            "outer_forward": outer_forward,
            "outer_reverse": outer_reverse,
            "inner_forward": inner_forward,
            "inner_reverse": inner_reverse
        },
        assay_folder
    )

    # -----------------------------
    # Write ApE feature library
    # -----------------------------
    feature_file = write_ape_feature_library(
        safe_id,
        outer_forward,
        outer_reverse,
        inner_forward,
        inner_reverse,
        assay_folder
    )

    return primer_file, csv_file, ape_file, feature_file


# -----------------------------
# Batch assay design
# -----------------------------
MAX_ATTEMPTS_PER_ASSAY = 10

def batch_assay_pipeline(
    organisms,
    assays_per_org,
    region_length,
    email,
    outer_primer_settings,
    inner_primer_settings
):

    # Use these settings directly instead of re-reading the config
    outer_settings = outer_primer_settings or {}
    inner_settings = inner_primer_settings or {}

    fetcher = NCBISequenceFetcher(email=email)

    for org_name, accession in organisms.items():

        print(f"\nProcessing organism: {org_name}")

        genome_path, safe_id = fetcher.fetch_genome(accession)

        org_folder = OUTPUT_FOLDER / f"{safe_id}_{org_name}"
        org_folder.mkdir(exist_ok=True)

        record = next(SeqIO.parse(genome_path, "fasta"))

        genome_seq = str(record.seq)
        genome_len = len(genome_seq)

        for i in range(1, assays_per_org + 1):

            attempt = 0

            while attempt < MAX_ATTEMPTS_PER_ASSAY:

                attempt += 1

                start = random.randint(0, genome_len - region_length)
                end = start + region_length

                assay_seq = genome_seq[start:end]

                assay_folder = org_folder / f"{safe_id}_Assay_{i:03d}"
                assay_folder.mkdir(exist_ok=True)

                info_file = assay_folder / "info.txt"

                with open(info_file, "w") as f:
                    f.write(f"Organism: {org_name}\n")
                    f.write(f"Genome accession: {accession}\n")
                    f.write(f"Genome location: {start}-{end}\n")
                    f.write(f"Region length: {region_length}\n")

                assay_safe_id = f"{safe_id}_{start}_{end}"

                try:

                    primer_file, csv_file, ape_file, feature_file = design_assay(
                        assay_seq,
                        assay_safe_id,
                        assay_folder,
                        outer_settings,
                        inner_settings
                    )

                except KeyError:

                    print(f"Attempt {attempt}: Inner primers not found, retrying...")

                    assay_folder.rmdir()

                    continue

                print(
                    f"Assay {i:03d} completed: {assay_folder.name} "
                    "(TXT, CSV, ApE, Feature Library generated)"
                )

                break

            else:

                print(
                    f"Warning: Could not generate assay {i:03d} "
                    f"after {MAX_ATTEMPTS_PER_ASSAY} attempts"
                )
