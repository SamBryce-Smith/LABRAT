#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import pyfaidx
from functools import reduce
from Bio.Seq import Seq
import sys


def _n_ids(gr, id_col, is_df=False):

    assert id_col in gr.columns

    if not is_df:
        return len(set(gr.as_df()[id_col]))
    else:
        return gr[id_col].nunique()


# get_terminal_regions (last n)
# add_region_number
    # add a '0-based' parameter (to count region_number from 0)


def _rev_complement_seq(df, seq_col="seq"):
    '''
    Reverse complement sequence in Seq objects if region found on minus strand
    Use internally inside a pr.assign/apply
    Returns a Series of Seq objects (unmodified if on plus strand, complemented if on minus)
    '''
    if (df["Strand"] == "+").all():

        return df[seq_col]

    elif (df["Strand"] == "-").all():

        return df[seq_col].apply(lambda seq: seq.reverse_complement())

    else:
        raise ValueError("Invalid value in 'Strand' column for df - must be one of '+' or '-'")


def get_terminal_regions(gr,
                         feature_col="Feature",
                         feature_key="exon",
                         id_col="transcript_id",
                         region_number_col="exon_number",
                         region_number_type="strand_specific",
                         ):
    '''
    Return gr of last exons for each transcript_id
    In process, region_number_col will be converted to type 'int'
    StringTie merged GTFs (or whatever tool single_steps/stringtie_longreads.smk is using)
    reports exon_number that DOES NOT RESPECT STRAND (from browsing in IGV)
    i.e. for minus-strand - largest exon_number for transcript corresponds to FIRST EXON, not last
    Annotated (i.e. Ensembl) reported exon_numbers DO RESPECT STRAND (i.e. max always = last exon)

    if Do respect strand, put source = None (default)
    if Don't respect strand, put source = "stringtie" (i.e. plus strand = max, minus strand = min)
    '''

    assert region_number_type in ["strand_specific", "not_strand_specific"]
    # assert which_region in ["first", "last"]
    assert region_number_col in gr.columns.tolist()
    assert feature_col in gr.columns.tolist()
    assert id_col in gr.columns.tolist()

    # Make sure only 'exon' features are in the gr
    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)


    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    mod_gr = gr.apply(lambda df: df.sort_values(by=[id_col, region_number_col], ascending=True),
                      nb_cpu=1)


    if region_number_type == "strand_specific":
        # strand_specific means that 1 = first region of group regardless of strand
        # Pick last region entry by max region number for each transcript (id_col)
        # Pick first region entry by min region number for each transcript (id_col)
        # keep="last" sets last in ID to 'False' and all others true (negate to keep last only)
        # keep="first" sets first in ID to 'False'

        out_gr = mod_gr.subset(lambda df: ~(df.duplicated(subset=[id_col], keep="last")),
                               nb_cpu=1
                               )

    elif region_number_type == "not_strand_specific":
        # Numbering Doesn't respect strand
        # Need to flip selecting first/last in group depending on strand
        # minus strand - pick min if Minus strand, max if plus strand

        # + strand - pick last in group/ID
        # - strand - pick first in group/ID
        out_gr = (mod_gr.subset(lambda df:
                                #1. plus strand & last in group/ID
                                (df["Strand"] == "+") & ~(df.duplicated(subset=[id_col],
                                                                        keep="last")) |
                                #2. minus strand & first in group/ID
                                (df["Strand"] == "-") & ~(df.duplicated(subset=[id_col],
                                                                        keep="first")),
                                nb_cpu=1)
                 )


    return out_gr



def get_last_two_exons(gr, id_col, region_id_col):
    '''
    '''

    # assert

    # First extract last exons for each ID
    last_exons = get_terminal_regions(gr, id_col=id_col, region_number_type="strand_specific")

    # Get these region IDs, remove them from input gr
    last_ids = set(last_exons.as_df()[region_id_col])

    nl_gr = gr.subset(lambda df: ~df[region_id_col].isin(last_ids))

    # Another round of extracting last exons from gr with last exons removed
    penul_exons = get_terminal_regions(nl_gr,
                                       id_col=id_col,
                                       region_number_type="strand_specific")

    last_two = pr.concat([last_exons, penul_exons])

    return last_two


def main(in_gtf, in_fasta, out_prefix):
    '''
    '''
    # If have multiple 'tag' key-value pairs, values are stored as comma-separated strings in a single column
    gtf = pr.read_gtf(in_gtf, duplicate_attr=True)

    # Only need gene_id, transcript_id, gene_type, transcript_type & tag for analysis
    # subset to these cols to save memory
    gtf = gtf[["Feature",
               "gene_id",
               "transcript_id",
               "gene_type",
               "transcript_type",
               "exon_number",
               "tag"]]

    # Subset for protein-coding genes
    print("subset for protein_coding genes")
    print(_n_ids(gtf, "transcript_id"))

    gtf = gtf.subset(lambda df: df["gene_type"] == "protein_coding", nb_cpu=1)

    print(_n_ids(gtf, "transcript_id"))

    # Subset for protein-coding transcripts
    print("subset for protein_coding transcripts")

    gtf = gtf.subset(lambda df: df["transcript_type"] == "protein_coding",
                     nb_cpu=1
                     )

    print(_n_ids(gtf, "transcript_id"))

    # Filter transcripts with 'tag' - mRNA_end_NF
    print("subset for not having 'mRNA_end_NF'")
    gtf = gtf.subset(lambda df: ~df["tag"].str.contains("mRNA_end_NF", regex=False))

    print(_n_ids(gtf, "transcript_id"))

    # Additional filters = tx has at least two exons, length of terminal fragment > 200
    # At this point only need exons
    exons = gtf.subset(lambda df: df["Feature"] == "exon")

    # tx has at least two exons

    # Sort by tx_id & exon number?

    # Setting keep 'False' marks all duplicates as True - keeps transcripts with > 1 exon)
    print("subset for at least two exons")

    exons = exons.subset(lambda df: df.duplicated(subset=["transcript_id"], keep=False), nb_cpu=1)

    print(_n_ids(exons, "transcript_id"))

    # extract last two exons for each transcript
    exons = exons.assign("exon_number",
                         lambda df: df["exon_number"].astype(float).astype(int),
                         nb_cpu=1)

    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    exons = exons.apply(lambda df: df.sort_values(by=["transcript_id", "exon_number"],
                                                  ascending=True),
                        nb_cpu=1)

    exons = exons.assign("region_id",
                         lambda df: df["transcript_id"].str.cat(df[["Start", "End", "Strand"]].astype(str),
                                                                sep=':')
                         )

    last_two_exons = get_last_two_exons(exons, id_col="transcript_id", region_id_col="region_id")
    # print(last_two_exons.transcript_id.value_counts().loc[lambda x: x != 2])
    # print(last_two_exons.apply(lambda df: df.sort_values(by=["transcript_id", "exon_number"],
                                          # ascending=True),
                                          # nb_cpu=1).head(n=20))

    # Filter for transcripts where terminal fragment size (length of last two exons) is > 200nt
    # lengths of each exon as a column
    last_two_exons.exon_length = last_two_exons.lengths()

    # iterable of set of tx_ids passing length filter for each chr/strand tuple
    term_lengths = (last_two_exons.apply(lambda df:
                                         (set(df.groupby("transcript_id")
                                              ["exon_length"].sum()
                                              # subset series with tx_id as index
                                              .loc[lambda x: x > 200]
                                              .index)
                                              ),
                                         as_pyranges=False
                                         ).values()
                    )

    # print(term_lengths)

    min_len_ids = reduce(lambda x, y: x.union(y), term_lengths)

    # Subset for min length transcripts
    print("subsetting for min length of fragment")
    last_two_exons = last_two_exons.subset(lambda df: df["transcript_id"].isin(min_len_ids))

    print(_n_ids(last_two_exons, "transcript_id"))

    # Get sequence for each terminal fragment
    # chr strand tuple invalid with std FASTAs
    # Do rev comp myself
    seq = pr.get_fasta(last_two_exons.unstrand(), in_fasta)
    seq = seq.apply(lambda x: Seq(x))
    last_two_exons.seq = seq

    last_two_exons = last_two_exons.assign("seq", lambda df: _rev_complement_seq(df))
    # Back to string objects
    last_two_exons = last_two_exons.assign("seq", lambda df: df["seq"].apply(lambda x: str(x)))


    # last_two_exons.seq = pr.get_fasta(last_two_exons, in_fasta)

    # Sequences are per exon, need to concatenate

    # Make sure gr is sorted by transcript_id & 'region number' (ascending order so 1..n)
    # When concat by tx exons are in correct order
    last_two_exons = last_two_exons.apply(lambda df: df.sort_values(by=["transcript_id", "exon_number"],
                                          ascending=True),
                                          nb_cpu=1)

    # Get dict of {(chr, strand): transcript_id | seq }
    tx_seqs = last_two_exons.apply(lambda df: (df.groupby("transcript_id")
                                               .apply(lambda df: df["seq"].str.cat())
                                               .reset_index()
                                               .rename(columns={0: "seq"})
                                               ),
                                    as_pyranges=False
                                   )

    # Combine into single_df
    tx_seqs = pd.concat(tx_seqs.values())

    #Write to FASTA
    with open(out_prefix + ".fa", "w") as outfa:
        for _, row in tx_seqs.iterrows():
            outfa.write(">" + row["transcript_id"].split(".")[0] + "\n" + row["seq"] + "\n")

    print("Fasta done!")

    last_two_exons.drop(["seq", "region_id"]).to_gtf(out_prefix + ".gtf")
    print("GTF done!")
    # exons.apply(lambda df: df.loc[df.groupby("transcript_id")["exon_number"].nlargest(2), :])


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
