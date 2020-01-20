import sys
def convertuniprobe(mot):
    unidict={}
    with open("/home/kipkurui/Project/Motif_Assessment/Data/uniprobemouse.txt") as uni:
        for line in uni:
            unidict[line.split()[1].lower()]=line.split()[0]
    return unidict[mot.lower()]
mot=sys.argv[1]
print convertuniprobe(mot)
