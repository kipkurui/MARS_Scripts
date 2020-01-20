"""
	libkmersvm.py; common library for kmersvm_train.py and kmersvm_classify.py
	Copyright (C) 2011 Dongwon Lee

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os
import os.path
import optparse

from bitarray import bitarray


def bitarray_fromfile(filename):
    """
    """
    fh = open(filename, 'rb')
    bits = bitarray()
    bits.fromfile(fh)

    return bits, fh


def generate_kmers(kmerlen):
    """make a full list of k-mers

    Arguments:
    kmerlen -- integer, length of k-mer

    Return:
    a list of the full set of k-mers
    """

    nts = ['A', 'C', 'G', 'T']
    kmers = []
    kmers.append('')
    l = 0
    while l < kmerlen:
        imers = []
        for imer in kmers:
            for nt in nts:
                imers.append(imer + nt)
        kmers = imers
        l += 1

    return kmers


def revcomp(seq):
    """get reverse complement DNA sequence

    Arguments:
    seq -- string, DNA sequence

    Return:
    the reverse complement sequence of the given sequence
    """
    rc = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
    return ''.join([rc[seq[i]] for i in xrange(len(seq) - 1, -1, -1)])


def generate_rcmap_table(kmerlen, kmers):
    """make a lookup table for reverse complement k-mer ids for speed

    Arguments:
    kmerlen -- integer, length of k-mer
    kmers -- list, a full set of k-mers generated by generate_kmers

    Return:
    a dictionary containing the mapping table
    """
    revcomp_func = revcomp

    kmer_id_dict = {}
    for i in xrange(len(kmers)):
        kmer_id_dict[kmers[i]] = i

    revcomp_mapping_table = []
    for kmerid in xrange(len(kmers)):
        rc_id = kmer_id_dict[revcomp_func(kmers[kmerid])]
        if rc_id < kmerid:
            revcomp_mapping_table.append(rc_id)
        else:
            revcomp_mapping_table.append(kmerid)

    return revcomp_mapping_table


def read_fastafile(filename, subs=True):
    """Read sequences from a file in FASTA format

    Arguments:
    filename -- string, the name of the sequence file in FASTA format
    subs -- bool, substitute 'N' with 'A' if set true

    Return:
    list of sequences, list of sequence ids
    """

    sids = []
    seqs = []

    try:
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()

    except IOError, (errno, strerror):
        print "I/O error(%d): %s" % (errno, strerror)
        sys.exit(0)

    seq = []
    for line in lines:
        if line[0] == '>':
            sids.append(line[1:].rstrip('\n').split()[0])
            if seq != []: seqs.append("".join(seq))
            seq = []
        else:
            if subs:
                seq.append(line.rstrip('\n').upper().replace('N', 'A'))
            else:
                seq.append(line.rstrip('\n').upper())

    if seq != []:
        seqs.append("".join(seq))

    return seqs, sids
