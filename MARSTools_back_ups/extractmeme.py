import sys

if len(sys.argv) < 3:
    print("\nUsage: python extractmeme.py <meme-file-input> <meme-out> <motif> \n")
    sys.exit(1)
meme = sys.argv[1]
# output=sys.argv[2]
motif = sys.argv[2]
# memehead=sys.argv[3]
# motifname=sys.argv[4]

with open(meme) as f1:
    # with open(output, "w") as f2:
    Lines = f1.readlines()
    print """MEME version 4.4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000 
""",
    for i, line in enumerate(Lines):
        head = line.split()

        if motif in line and motif == head[1]:
            k = i
            print ("\n" + Lines[i])
            if "log-odds" in Lines[i + 1]:
                odds = Lines[k + 2].split()
                # print odds[5]
                for j in range(2, (int(odds[5]) + 3) * 2):
                    # f2.write(Lines[i+j])
                    print (Lines[i + j]),
            elif "letter-probability" in Lines[i + 2]:
                print Lines[i + 2],
                odds = Lines[k + 2].split()
                # print odds
                # print Lines[i+1],
                # f2.write(odds[5])
                for j in range(0, (int(odds[5]))):
                    # f2.write(Lines[i+3+j])
                    print (Lines[i + 3 + j]),
            else:
                print"No data"
                sys.exit(1)
