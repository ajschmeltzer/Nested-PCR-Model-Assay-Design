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
base_folder = Path(__file__).parent.parent.resolve()
genome_folder = base_folder / "Designed Assays" / "Genomes"
output_folder = base_folder / "Designed Assays" / "Assays by Organism"

# Make sure directories exist
genome_folder.mkdir(parents=True, exist_ok=True)
output_folder.mkdir(parents=True, exist_ok=True)

#Set base folder
def get_base_folder():
    return Path(os.environ.get("PCR_ASSAY_ROOT", base_folder))

# -----------------------------
# Check for existing assay folders
# -----------------------------
# Finds the next available assay number by scanning existing assay folders
# (e.g., <safe_id>_Assay_001) and continuing the sequence

def get_next_assay_number(org_folder, safe_id):
    org_folder.mkdir(parents=True, exist_ok=True)
    existing_numbers = []
    pattern = re.compile(rf"{re.escape(safe_id)}_Assay_(\d+)") #Folder name pattern for the organism
    
    #Checking for existing assay folders
    for f in org_folder.iterdir():
        if f.is_dir():
            match = pattern.fullmatch(f.name)
            if match:
                existing_numbers.append(int(match.group(1)))
    
    #Goes to next available number            
    if existing_numbers:
        next_number = max(existing_numbers) + 1
    else:
        next_number = 1
        
    return next_number

# -----------------------------
# Helper functions
# -----------------------------
#Calculates GC content of a sequence
def gc_content(seq):
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return 100 * gc_count / len(seq)

#Creates the reverse compliment of a sequence
def revcomp(seq):
    complement = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(complement)[::-1]

# -----------------------------
# NCBI fetcher
# -----------------------------
# Downloads the genome region from the config file
# Skips the download if the genone region is already present in the 'Genomes' folder

class NCBISequenceFetcher:
    def __init__(self, email, genome_dir=None):
        Entrez.email = email #entrez requires a user email to access the database
        #Making folder to store specific genomes in
        if genome_dir is None:
            genome_dir = genome_folder
        self.genome_dir = Path(genome_dir)
        self.genome_dir.mkdir(parents=True, exist_ok=True)

    def fetch_genome(self, accession):
        genome_path = self.genome_dir / f"{accession}.fasta"
        if genome_path.exists():
            print(f"Genome already downloaded: {genome_path.name}")
            safe_id = re.sub(r'[:<>"/\\|?*]', '_', accession) #Need to remove the ":" from the name of the NCIB genome range
            return genome_path, safe_id

        print(f"Fetching full genome from NCBI: {accession} ...")
        handle = Entrez.efetch(
            db="nuccore",
            id=accession,
            rettype="fasta",
            retmode="text",
            #mask="hard"  # Downloads hard masked 'N' genome. This doesn't actually work. Will need a different solution for masking
        )
        record = SeqIO.read(handle, "fasta")
        SeqIO.write(record, genome_path, "fasta")
        safe_id = re.sub(r'[:<>"/\\|?*]', '_', accession)
        print(f"Genome saved: {genome_path}")
        return genome_path, safe_id

# -----------------------------
# ApE file writer
# -----------------------------
# Function for writing the outer amplicon sequence to an ApE file
# Annotates the ApE file with the designed primers

def write_ape_file(assay_name, outer_amp_seq, primers, out_dir):

    #Generic placeholder primer names
    primer_colors = {
        "outer_forward": "#FF0000",   # red
        "outer_reverse": "#FF0000",
        "inner_forward": "#00FFFF",   # cyan
        "inner_reverse": "#00FFFF"
    }

    #Naming the ApE file
    ape_file = out_dir / f"{assay_name}.ape"
    seq_upper = outer_amp_seq.upper() #uppercase letters

    with open(ape_file, "w") as f:

        f.write(f"LOCUS       {assay_name} {len(outer_amp_seq)} bp DNA linear\n")
        f.write(f"DEFINITION  Outer amplicon for {assay_name}\n")
        f.write("FEATURES             Location/Qualifiers\n")

        for name, primer_seq in primers.items():

            primer_seq = primer_seq.upper()
            rc = revcomp(primer_seq)
            color = primer_colors[name]

            pos = seq_upper.find(primer_seq)

            #Annotating the ApE file with the designed primers
            if pos != -1:
                start = pos + 1
                end = pos + len(primer_seq)
                f.write(f"     primer_bind     {start}..{end}\n")
                f.write(f"                     /label=\"{name}\"\n")
                f.write(f"                     /ApEinfo_fwdcolor={color}\n")
                f.write(f"                     /ApEinfo_revcolor={color}\n")

            #Reverse Primers
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

        # Writing the outer amplicon sequence in the ApE file
        # Formatting sequence in chunks of 60 bases, needed for ApE formatting
        for i in range(0, len(outer_amp_seq), 60):
            seq_chunk = outer_amp_seq[i:i+60].upper()
            f.write(f"{i+1:9} {seq_chunk}\n")

        f.write("//\n")

    return ape_file

# -----------------------------
# ApE feature library writer
# -----------------------------
# Function for writing a feature library in the ApE format for future use

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
# Designs nested PCR primers for a given sequence using settings from the config file
# Saves the primer information in a .csv and .txt file
# Writes the annotated ApE file and ApE feature library for the designed primers

def design_assay(sequence, safe_id, assay_folder, outer_settings, inner_settings):
    
    # -----------------------------
    # Outer primers
    # -----------------------------
    #Designing outer primers
    outer_results = primer3.bindings.design_primers(
        {'SEQUENCE_ID': safe_id, 'SEQUENCE_TEMPLATE': sequence},
        outer_settings
    )

    #Outer primer results
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
    outer_amplicon_tm = primer3.calc_tm(outer_amplicon_seq) #Amplicon Tm based on GC content, not very accurate

    # -----------------------------
    # Inner primers
    # -----------------------------
    # Define inner assay region between outer primers with a 5 bp buffer on each side
    inner_start = outer_left_start + outer_left_len + 5
    inner_end = outer_right_start - outer_right_len - 5
    inner_sequence = sequence[inner_start:inner_end]

    #Designing inner primers
    inner_results = primer3.bindings.design_primers(
        {'SEQUENCE_ID': safe_id + "_inner", 'SEQUENCE_TEMPLATE': inner_sequence},
        inner_settings
    )

    #Inner primer results
    inner_forward = inner_results['PRIMER_LEFT_0_SEQUENCE']
    inner_reverse = inner_results['PRIMER_RIGHT_0_SEQUENCE']
    inner_product_size = inner_results['PRIMER_PAIR_0_PRODUCT_SIZE']
    inner_left_start, inner_left_len = inner_results['PRIMER_LEFT_0']
    inner_amplicon_seq = inner_sequence[inner_left_start : inner_left_start + inner_product_size]

    inner_forward_tm = primer3.calc_tm(inner_forward)
    inner_reverse_tm = primer3.calc_tm(inner_reverse)
    inner_forward_gc = gc_content(inner_forward)
    inner_reverse_gc = gc_content(inner_reverse)
    inner_amplicon_tm = primer3.calc_tm(inner_amplicon_seq) #Amplicon Tm based on GC content, not very accurate

    # ----------------------------- 
    # Save TXT
    # -----------------------------
    primer_file = assay_folder / f"{safe_id}_primers.txt"

    with open(primer_file, "w") as f:
        f.write("=== Outer Primers ===\n")
        f.write(f"Forward: {outer_forward.upper()}\n")
        f.write(f"Reverse: {outer_reverse.upper()}\n")
        f.write(f"Reverse (5'->3' revcomp): {revcomp(outer_reverse).upper()}\n")
        f.write(f"Forward Tm: {outer_forward_tm:.2f} °C\n")
        f.write(f"Reverse Tm: {outer_reverse_tm:.2f} °C\n")
        f.write(f"Forward GC%: {outer_forward_gc:.2f}\n")
        f.write(f"Reverse GC%: {outer_reverse_gc:.2f}\n")
        f.write(f"Product size: {outer_product_size}\n")
        f.write(f"Amplicon Tm (GC approx.): {outer_amplicon_tm:.2f} °C\n")
        f.write(f"Amplicon sequence: {outer_amplicon_seq.upper()}\n\n")

        f.write("=== Inner Primers ===\n")
        f.write(f"Forward: {inner_forward.upper()}\n")
        f.write(f"Reverse: {inner_reverse.upper()}\n")
        f.write(f"Reverse (5'->3' revcomp): {revcomp(inner_reverse).upper()}\n")
        f.write(f"Forward Tm: {inner_forward_tm:.2f} °C\n")
        f.write(f"Reverse Tm: {inner_reverse_tm:.2f} °C\n")
        f.write(f"Forward GC%: {inner_forward_gc:.2f}\n")
        f.write(f"Reverse GC%: {inner_reverse_gc:.2f}\n")
        f.write(f"Product size: {inner_product_size}\n")
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
            "Outer", 
            outer_forward.upper(), 
            outer_reverse.upper(), 
            revcomp(outer_reverse).upper(),
            
            f"{outer_forward_tm:.2f}", 
            f"{outer_reverse_tm:.2f}",
            f"{outer_forward_gc:.2f}", 
            f"{outer_reverse_gc:.2f}",
            outer_product_size,
            f"{outer_amplicon_tm:.2f}", 
            outer_amplicon_seq.upper()
        ])

        writer.writerow([
            "Inner", 
            inner_forward.upper(), 
            inner_reverse.upper(), 
            revcomp(inner_reverse).upper(),
            f"{inner_forward_tm:.2f}", 
            f"{inner_reverse_tm:.2f}",
            f"{inner_forward_gc:.2f}", 
            f"{inner_reverse_gc:.2f}",
            inner_product_size,
            f"{inner_amplicon_tm:.2f}", 
            inner_amplicon_seq.upper()
        ])

    # -----------------------------
    # Write annotated ApE file
    # -----------------------------
    ape_file = write_ape_file(
        safe_id,
        outer_amplicon_seq,
        {
            #Placeholder generic primer names
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
#Calls all previously defined functions to design primers for the organisms provided in the config file

#Sometimes Primer3 is unable to design inner with the given config settings.
# If this occurs, a new region of the genome will be designed for up to the max attempts value
max_attempts_per_assay = 10

#Calling user settings from the config file
def batch_assay_pipeline(
    organisms,
    assays_per_org,
    region_length,
    email,
    outer_primer_settings,
    inner_primer_settings
):

    # Settings from the config file
    outer_settings = outer_primer_settings or {}
    inner_settings = inner_primer_settings or {}

    #Calls the genome fetcher class
    fetcher = NCBISequenceFetcher(email=email)

    #Iterates through the organisms in the config file
    for org_name, accession in organisms.items():

        print(f"\nProcessing organism: {org_name}")

        #Fetching genomes and storing them in the genome folder
        genome_path, safe_id = fetcher.fetch_genome(accession)

        #Creating specific organsim folder
        org_folder = output_folder / f"{safe_id}_{org_name}"
        org_folder.mkdir(exist_ok=True)

        #Opening the fasta sequence file
        record = next(SeqIO.parse(genome_path, "fasta"))

        genome_seq = str(record.seq)
        genome_len = len(genome_seq)

        #Creates specified number of assays for each organism
        next_i = get_next_assay_number(org_folder, safe_id)
        for i in range(next_i, next_i + assays_per_org):

            attempt = 0

            while attempt < max_attempts_per_assay:

                attempt += 1

                #Assays are designed for random regions of the genome. A random number generator is used to pick the regions
                start = random.randint(0, genome_len - region_length) 
                end = start + region_length

                assay_seq = genome_seq[start:end]

                #Creating specific assay folder
                assay_folder = org_folder / f"{safe_id}_Assay_{i:03d}"
                assay_folder.mkdir(exist_ok=True)

                #Info file of the genome location is recorded
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

                #If an assay design fails
                except KeyError:
                    print(f"Attempt {attempt}: Inner primers not found, retrying...")

                    #Removing directory for failed assay design
                    assay_folder.rmdir()

                    continue

                print(
                    f"Assay {i:03d} completed: {assay_folder.name} "
                    "(TXT, CSV, ApE, Feature Library generated)"
                )

                break
            
            #If no assays can be designed within the max number of attempts
            else:
                print(
                    f"Warning: Could not generate assay {i:03d} "
                    f"after {max_attempts_per_assay} attempts"
                )
