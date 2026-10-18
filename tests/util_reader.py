def fasta_names(path):
    with open(path) as fh:
        return {line[1:].split()[0] for line in fh if line.startswith(">")}
