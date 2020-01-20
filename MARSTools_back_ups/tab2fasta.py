import sys


def tab2fasta(posneg, fasta):
    i = 0
    # print fasta
    posneglen = sum(1 for line in open(posneg))
    # close(posneg)
    with open(posneg) as tab:
        with open(fasta, 'w') as fa:
            for line in tab:
                details = line.split()
                if len(details) == 2:
                    pos = 1
                else:
                    pos = 2
                if i < posneglen / 2:  # This has to be changed to consider the sequences with sequences less than 500
                    fa.write(">" + line.split()[0] + '\n' + line.split()[pos] + "\n")
                i += 1


tab2fasta(sys.argv[1], sys.argv[2])
