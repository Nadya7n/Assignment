import argparse
from Bio import SeqIO
from tqdm import tqdm


def main(settings):

    ref_path = settings["fasta"]
    target_path = settings["target"]
    output_path = settings["output"]
    counter = 0

    with open(target_path) as fh:
        with open(ref_path) as fh_2:
            genome = [str(record.seq).lower() for record in SeqIO.parse(fh_2, "fasta")]
            for idr, record_2 in enumerate(tqdm(SeqIO.parse(fh, "fasta"))):
                for chromosome in genome:
                    target_sequence = str(record_2.seq).lower()
                    counter += chromosome.count(target_sequence)
                if counter > 1:
                    find_homo(ref_path, output_path, target_sequence, record_2.id)
                counter = 0


def find_homo(ref, output_path, target_seq, target_id):
    with open(output_path, "a+") as fw:
        target_seq = target_seq.lower()

        first_sep = target_id.find(":")
        second_sep = target_id.find("-")

        chromosome_target = target_id[:first_sep]
        start_target = int(target_id[first_sep + 1 : second_sep])
        end_target = int(target_id[second_sep + 1 :])
        length_target = end_target - start_target

        fw.write(f"{chromosome_target} {start_target} {end_target}\n")
        with open(ref) as fh:
            genome = [str(record.seq).lower() for record in SeqIO.parse(fh, "fasta")]
            for chrom in genome:
                first_match = chrom.find(target_seq, 0, start_target)
                second_match = chrom.find(target_seq, end_target)
                chromosome_ref = genome.index(chrom) + 1
                if first_match != -1 or second_match != -1:
                    start_ref = first_match if first_match != -1 else second_match
                    end_ref = start_ref + length_target
                    fw.write(f"chr{chromosome_ref} {start_ref} {end_ref}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Takes 3 arguments")

    parser.add_argument(
        "-f", "--fasta", help="Fasta file of reference genome", required=True
    )
    parser.add_argument(
        "-t", "--target", help="Fasta file of target sequences", required=True
    )
    parser.add_argument(
        "-o", "--output", help="Coordinates homologous sequences", required=True
    )

    args = vars(parser.parse_args())

    settings = {
        "fasta": args["fasta"],
        "target": args["target"],
        "output": args["output"],
    }

    main(settings)
