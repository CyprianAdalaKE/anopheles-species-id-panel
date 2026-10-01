import os
import sys

# example:
# >AAAB01008983|kraken:taxid|7165 | organism=Anopheles_gambiae_PEST | version=AgamP4 | length=28559 | SO=supercontig |
# >AgamP4_2L|kraken:taxid|{}
# organism=Anopheles_gambiae_PEST
# version=AgamP4|length=49364325|SO=chromosome

def main(fasta_file: str, taxon_id: int, output_dir: str) -> None:
    """Add a taxonmy ID to a file"""

    new_fasta_file = output_dir + "/" + os.path.basename(fasta_file).replace(".fasta", ".k2.fasta")
    with open(new_fasta_file, "w") as new_fasta:
        with open(fasta_file, "r") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    fields = line.replace(" ", "").split("|")
                    new_header = fields[0] + "|kraken:taxid|" + taxon_id + "\n"
                    print(new_header)
                    new_fasta.write(new_header)
                else:
                    new_fasta.write(line)


if __name__ == "__main__":
    main(fasta_file=sys.argv[1], 
         taxon_id=sys.argv[2],
         output_dir=sys.argv[3])


