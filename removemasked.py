"""
removemasked.py removes the sequences that have been repeat masked

Takes as imput fasta sequence with repeat masks and writes out a fasta
withot repeat masked sequences.

Usage:
    python removemasked.py <input-fasta> <output-fasta>

"""
import sys


def removemasked(fa, out):
    '''
    Removes fasta sequences that have been repeatmasked
    '''
    with open(fa) as fas:
        with open(out, "w") as faout:
            for line in fas:
                if "N" in line:
                    continue
                else:
                    faout.write(line)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print __doc__
        sys.exit(1)
    fa = sys.argv[1]
    out = sys.argv[2]
    removemasked(fa, out)
