from Bio import SeqIO
import argparse


def main(settings):
    """
    Change headers in fasta file
    :return: suitable headers
    """

    input_path = settings["input"]
    output_path = settings["output"]
    target_word = "chromosome"

    with open(input_path) as fh:
        with open(output_path, "a+") as fw:
            for record in SeqIO.parse(fh, "fasta"):
                description = record.description.split()
                if target_word in description and description[0].startswith("NC"):
                    index_of_chr_number = description.index(target_word) + 1
                    number_chr = description[index_of_chr_number].strip(",")
                    fw.write(f">chr{number_chr}\n{record.seq}\n")


if __name__ == "__main__":
    """
    Takes arguments from the command line,
    runs the main function with the accepted arguments
    """
    
    parser = argparse.ArgumentParser(description="Takes 2 arguments")

    parser.add_argument("-i", "--input", help="Fasta input file", required=True)
    parser.add_argument("-o", "--output", help="Fasta output file", required=True)

    args = vars(parser.parse_args())

    settings = {"input": args["input"], "output": args["output"]}

    main(settings)
