import re

IUPAC_CODES = {
    'A':    'A',
    'C':    'C',
    'G':    'G',
    'T':    'T',
    'U':    'U',
    'R':    '[AG]',
    'Y':    '[CTU]',
    'S':    '[GC]',
    'W':    '[ATU]',
    'K':    '[GTU]',
    'M':    '[AC]',
    'B':    '[CGTU]',
    'D':    '[AGTU]',
    'H':    '[ACTU]',
    'V':    '[ACG]',
    'N':    '[ACGTU]',
}

def load_fasta(filename):
    sequences = {}

    current = None
    for line in open(filename):
        if line.startswith('>'):
            current = line[1:].split(None, 1)[0]
            sequences[current] = ''
        else:
            sequences[current] += line.strip()

    return sequences

def iupac_to_regexp(iupac):
    return re.compile(''.join(IUPAC_CODES[code] for code in iupac))

def load_restriction_enzyme_list(filename):
    resites = {}

    for line in open(filename):
        name, pattern = line[:-1].split(None, 1)
        resites[name] = iupac_to_regexp(pattern)

    return resites

if __name__ == '__main__':
    #print(load_fasta('human-miRNA-hairpins.fa'))
    print(load_restriction_enzyme_list('restriction-enzymes.txt'))

