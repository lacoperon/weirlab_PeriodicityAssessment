import re
import csv

def read_fasta_and_filter(fasta_file, gene_file):
    with open(fasta_file) as f:
        # Opens the gene_file, corresponding to all genes that have widespread
        # alterations in mRNA translation within 60 minutes
        gene_symbols = read_gene_file(gene_file)

        fasta_df = []
        header = ""
        sequence = ""
        hit_counter = 0
        gene_symbol = ""

        # This variable concatenates the dataset of sequences, obtaining
        # background nucleotide frequencies for the human genome
        total_sequence_aggregator = ""

        # Initializes dict containing counts for each gene_symbol
        gene_symbol_counter = {}
        for gene in gene_symbols:
            gene_symbol_counter[gene] = 0

        # Reads in all of lines in FASTA file
        while True:
            line = f.readline()
            # End condition
            if not line:
                # if we're exiting the loop, there should be an unwritten
                # header, sequence pair left to write
                if header and sequence and gene_symbol.strip() in gene_symbols:
                    # We write the unwritten pair, if it should be written
                    fasta_df.append([header, sequence, gene_symbol])
                    hit_counter += 1
                # Then we exit the loop
                break

            # Encountering a new FASTA header
            if line[0] == ">":
                # If appropriate, we write the header and sequence to our list
                if header and sequence and gene_symbol.strip() in gene_symbols:
                    fasta_df.append([header, sequence, gene_symbol])
                    hit_counter += 1
                    gene_symbol_counter[gene_symbol] += 1

                # We read in the new header
                header = line[1:]
                # Parsing the gene_symbol using RegEx
                gene_symbol = re.search("(?<=gene_symbol:)[^\s]*", header).group(0)
                # Initialize new empty sequence
                sequence = ""
            # Adding to sequence
            else:
                sequence += re.sub("\s", "", line)
                total_sequence_aggregator += re.sub("\s", "", line)

        # Printing out metadata for run to stdout
        print("Number of hits is {}".format(hit_counter))
        print("Average is {}".format(sum(gene_symbol_counter.values())/len(gene_symbol_counter.values())))
        print("Max is {}".format(max(gene_symbol_counter.values())))
        print("Min is {}".format(min(gene_symbol_counter.values())))
        with open("./data/hsapiens/runs/hits.csv", "w") as outfile:
           writer = csv.writer(outfile)
           writer.writerow(["num_hits"])
           writer.writerows(map(lambda x : [str(x)], gene_symbol_counter.values()))
           print("Wrote hits to csv, so we can visualize later")

        return [fasta_df, total_sequence_aggregator]

def read_gene_file(gene_file):
    with open(gene_file) as f:
        gene_names = []
        f.readline() # throws away header line
        while True:
            line = f.readline()
            if line:
                gene_names.append(line.strip())
            else:
                break

        return gene_names

def write_fasta_file(fasta_df):
    with open("./data/hsapiens/Homo_sapiens.cds.andreev.fa", "w") as outfile:
        writer = csv.writer(outfile)
        for fasta in fasta_df:
            writer.writerow(fasta)
    return

def calculate_nucleotide_freqs(total_seq):
    count_dict = {}
    for letter in total_seq:
        if letter not in count_dict:
            count_dict[letter] = 0
        count_dict[letter] += 1

    del count_dict["N"]

    total = count_dict["A"] + count_dict["T"] + count_dict["C"] + count_dict["G"]
    freq_dict = {k: v / float(total) for k, v in count_dict.items()}
    print(freq_dict.items())

fasta_file = "./data/hsapiens/Homo_sapiens.GRCh38.cds.all.fa"
gene_file  = "./data/hsapiens/andreev_genes.csv"
print(" --- Reading FASTA and filtering only for D.E. Andreev genes ---")
[fasta_df, total_seq] = read_fasta_and_filter(fasta_file, gene_file)
print(" --- Writing csv with FASTA values for those filtered genes  --- ")
write_fasta_file(fasta_df)
print(" --- Calculating background nucleotide frequencies for all genes ---")
calculate_nucleotide_freqs(total_seq)
