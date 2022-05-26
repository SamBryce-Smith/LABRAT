#!/usr/bin/env python3


from Bio import SeqIO
import sys

def main(me, them):

    me = SeqIO.to_dict(SeqIO.parse(me, 'fasta'))
    them = SeqIO.to_dict(SeqIO.parse(them, 'fasta'))

    print("Checking if transcripts in each are identical")
    print(sorted(list(me.keys())) == sorted(list(them.keys())))
    print("checking if sequences for each transcript are identical")
    print("tx_id added if sequences don't match")
    nm = [tx for tx in me.keys() if me[tx].seq != them[tx].seq]
    print(nm)
    # print(me)
    # print(me.keys())
    # print(str(me["ENST00000519638"].seq) == str(them["ENST00000519638"].seq))
    # print(me["ENST00000519638"].seq == them["ENST00000519638"].seq)
    # print(me["ENST00000519638"].seq[:50])
    # print(them["ENST00000519638"].seq[:50])



if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
