import mylib

HAIRPIN_FILE = 'human-miRNA-hairpins.fa'
RESENZ_FILE = 'restriction-enzymes.txt'

hairpin_seqs = mylib.load_fasta(HAIRPIN_FILE)
re_patterns = mylib.load_restriction_enzyme_list(RESENZ_FILE)

for name, pattern in re_patterns.items():
    print(f"Testing enzyme: {name}")

    matching_hairpins = [
        seqname
        for seqname, seq in hairpin_seqs.items()
        if pattern.search(seq)]

    print("  No of matching hairpins:", len(matching_hairpins))
    if len(matching_hairpins):
        print("  First few examples:", ', '.join(matching_hairpins[:5]))

    print()

