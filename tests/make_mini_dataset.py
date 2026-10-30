"""Usage: python3 make_mini_dataset.py OUTPUT_DIR

Write a tiny NanoASV input: a reference file, two barcode directories
and a metadata.csv. Each read is an exact copy of a reference sequence
(or of a sequence absent from the reference), so that every read has a
known fate. Sequences are random but seeded: the dataset is the same at
each call.
"""

import gzip
import os
import random
import sys

SEQUENCE_LENGTH = 1500  # within chopper's default length range
SEED = 42

rng = random.Random(SEED)


def random_sequence():
    return "".join(rng.choice("ACGT") for _ in range(SEQUENCE_LENGTH))


def reverse_complement(sequence):
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


# AB1.1.15 is a prefix of AB1.1.1500 (taxonomy lookups with grep -f)
REFERENCES = {
    "AB1.1.15": "Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus sp.",
    "AB1.1.1500": "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli",
    "REV1.1.1500": "Bacteria;Actinobacteriota;Actinobacteria;Micrococcales;Microbacteriaceae;Agromyces;Agromyces sp.",
    "MITO.1.1500": "Bacteria;Proteobacteria;Alphaproteobacteria;Rickettsiales;Mitochondria;Mitochondria;uncultured",
    "CHLO.1.1500": "Bacteria;Cyanobacteria;Cyanobacteriia;Chloroplast;Chloroplast;Chloroplast;uncultured",
    "EUK.1.1500": "Eukaryota;Amorphea;Obazoa;Opisthokonta;Holozoa;Choanozoa;uncultured",
    "SHORT.1.1500": "Bacteria;Planctomycetota;Phycisphaerae",
}
UNKNOWN = None

# per sample: (reference, number of reads, reverse strand)
SAMPLES = {
    "barcode01": [
        ("AB1.1.15", 3, False),
        ("REV1.1.1500", 2, True),
        ("MITO.1.1500", 1, False),
        ("CHLO.1.1500", 1, False),
        ("EUK.1.1500", 1, False),
        ("SHORT.1.1500", 1, False),
    ],
    "barcode02": [
        ("AB1.1.1500", 4, False),
        (UNKNOWN, 5, False),
    ],
}


def main():
    output_dir = sys.argv[1]
    sequences = {name: random_sequence() for name in REFERENCES}
    sequences[UNKNOWN] = random_sequence()

    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, "reference.fasta"), "w") as fasta:
        for name, taxonomy in REFERENCES.items():
            fasta.write(f">{name} {taxonomy}\n{sequences[name]}\n")

    for sample, reads in SAMPLES.items():
        os.makedirs(os.path.join(output_dir, sample), exist_ok=True)
        fastq_file = os.path.join(output_dir, sample, "reads.fastq.gz")
        with gzip.open(fastq_file, "wt") as fastq:
            read_number = 0
            for name, count, reverse in reads:
                sequence = sequences[name]
                if reverse:
                    sequence = reverse_complement(sequence)
                for _ in range(count):
                    read_number += 1
                    fastq.write(f"@{sample}_r{read_number}\n{sequence}\n"
                                f"+\n{'I' * SEQUENCE_LENGTH}\n")

    with open(os.path.join(output_dir, "metadata.csv"), "w") as metadata:
        metadata.write(",ID,Origin\n")
        for number, sample in enumerate(SAMPLES, start=1):
            metadata.write(f"{sample},{number},Mini\n")


if __name__ == "__main__":
    main()
