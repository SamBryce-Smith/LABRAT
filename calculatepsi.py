#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
import sys
from makeTFfasta import get_terminal_regions


def _df_add_region_number(df,id_col,sort_col="Start"):
    '''
    Return a Series of strand-aware region numbers (5'-3' in 1..n)
    Function to be used internally in a pr.assign (mainly by add_region_number)
    '''
    if id_col not in df.columns.tolist():
        raise KeyError(f"id_col - {id_col} - is not present in df for chr/strand pair {','.join([df.Chromosome.iloc[0], df.Strand.iloc[0]])}")

    elif (df.Strand == "+").all():
        # Start position smallest to largest = 5'-3'

        return df.groupby(id_col)[sort_col].rank(method="min", ascending=True)

    elif (df.Strand == "-").all():
        # Start position largest to smallest = 5'-3'

        return df.groupby(id_col)[sort_col].rank(method="min", ascending=False)

    elif df.empty:
        print("df is empty - returning empty pd.Series")
        return pd.Series()


def _pd_merge_gr(df, df_to_merge, how, on, suffixes, to_merge_cols):
    '''
    Perform a pd.merge inside a pr.apply to add columns from gr_to_merge based on metadata in gr_df
    Here, will check the chromosome and strand of provided gr (gr_df)
    and subset gr_to_merge to chr/strand before converting to df
    This should cut down on memory requirements when convert to df (which requires pd.concat() x chrom & strand)
    For this kind of merge, only expect joins between regions on same chr/strand
    '''
    #chromsomes returns a list of chr names (always 1 val)
    assert isinstance(to_merge_cols, list)

    if df_to_merge.empty:
        eprint("df_to_merge for chr/strand pair {} is empty - returning to_merge cols filled with NaNs".format(",".join([df.Chromosome.iloc[0], df.Strand.iloc[0]])))

        df_cols = df.columns.tolist()
        # on won't have suffix added - need to remove as a target column
        to_merge_cols = [col for col in to_merge_cols if col != on]

        # eprint(df_cols)
        # eprint(to_merge_cols)

        # list of cols shared between dfs - need to add suffixes[1]
        # list of cols only in df_to_merge - these stay the same in a merge
        only_cols = [col for col in to_merge_cols if col not in df_cols]
        shared_cols = [col for col in to_merge_cols if col in df_cols]

        # eprint("only_cols - {}".format(only_cols))
        # eprint("shared_cols - {}".format(shared_cols))

        out_shared = [col + suffixes[1] for col in shared_cols]
        target_cols = out_shared + only_cols

        nrows = len(df.index)

        out_cols = {col: pd.Series([np.nan]*nrows) for col in target_cols}

        return df.assign(**out_cols)

    else:
        return df.merge(df_to_merge,
                        how=how,
                        on=on,
                        suffixes=suffixes)


def cluster_to_posfactor(gr,
                         group_id_col="gene_id",
                         cluster_col="Cluster",
                         out_col="posfactor",
                         zero_based=True):
    '''
    Returns gr with 'out_col' column added
    where out_col is leftmost to rightmost cluster_col converted to a
    strand-aware 1..n order by group_id_col
    1 = most 5' site in group_id_col
    n = most 3' site in group_id_col
    '''

    # For groupby.rank to work appropriately, need a single row per 'Cluster'
    c2p = (gr[[group_id_col, cluster_col]]
           .apply(lambda df: df.drop_duplicates(subset=cluster_col)
                  )
           )

    # Add 1..n 5'-3' region number as a column
    c2p = c2p.assign(out_col,
                     lambda df: _df_add_region_number(df, group_id_col, cluster_col))

    # Return 'out_col' to original gr
    c2p_cols = c2p.columns.tolist()
    out_gr = gr.apply_pair(c2p, lambda df, df_to_merge: _pd_merge_gr(df,
                                                                     df_to_merge,how="left",
                                                                     on=cluster_col,
                                                                     suffixes=[None, "_match"],
                                                                     to_merge_cols=c2p_cols)
                           )

    # avoid adding extra 'PyRanges' cols (Chromosome etc.) from c2p
    out_gr = out_gr.drop(like="_match$")

    if zero_based:
        return out_gr.assign(out_col, lambda df: df[out_col] - 1)

    else:
        return out_gr



def add_region_number(gr,
                      id_col="transcript_id",
                      out_col="intron_number",
                      feature_key="intron",
                      feature_col="Feature",
                      sort_col="Start",
                      zero_based=False,
                      nb_cpu=1):
    '''
    Adds column to gr containing a strand aware region number column,
    ordered 5'-3' 1..n by a group of features (e.g. transcript)
    '''

    # Make sure only 'feature_key' rows are in the gr
    assert gr.as_df()[feature_col].drop_duplicates().tolist() == [feature_key], "only {} entries should be present in gr".format(feature_key)

    # Make sure sorted by position first.
    gr = gr.sort()

    # Add in region number column in strand aware manner, so 1 = most 5', n = most 3'

    gr = gr.assign(out_col, lambda df: _df_add_region_number(df, id_col, sort_col=sort_col), nb_cpu=nb_cpu)

    if not zero_based:
        return gr

    else:
        return gr.assign(out_col, lambda df: df[out_col] - 1)


def getpositionfactors(gr, cluster_distance, out_col="posfactor"):
    '''
    Group 3'ends of regions within merge window and assign 5'-3' position in gene
    '''

    end_3p = gr.three_end()

    # group 3'ends occurring within 25nt together with common ID
    end_3p = end_3p.cluster(by="gene_id",
                            slack=cluster_distance,
                            strand=True)

    # pr.cluster assigns an integer leftmost-rightmost per gene
    # i.e. is not strand specific - need to consider leftmost as 1st on plus ( & opposite for rightmost)
    # i.e. the smallest int per gene_id = most 5' in group
    end_3p = cluster_to_posfactor(end_3p,
                                  group_id_col="gene_id",
                                  out_col=out_col)

    return end_3p





def main(in_gtf, cluster_distance, out_prefix):
    '''
    '''

    # This should be output of makeTFfasta.py
    gtf = pr.read_gtf(in_gtf)

    gtf = gtf.assign("exon_number",
                     lambda df: df["exon_number"].astype(float).astype(int))

    # Get last exons of each tx (so can extract PAS)
    le = get_terminal_regions(le, id_col="transcript_id", region_number_type="strand_specific")

    # posfactor col added - zero-based 5'-3' position within gene
    pas = getpositionfactors(le, cluster_distance)

    # Can report positions of PAS/clusters later...
    merged_pas = pas.merge(slack=cluster_distance, by="gene_id")

    # Number of pas per gene (as dict of dfs)
    n_per_gene = pas.apply(lambda df: (df.groupby("gene_id")
                                       ["posfactor"]
                                       .nunique()
                                       .reset_index()
                                       .rename({"posfactor": "numberofposfactors"})),
                           as_pyranges=False
                           )

    # merge back n_per_gene into pas
    n_per_gene = pd.concat(n_per_gene.values())
    per_gene_cols = n_per_gene.columns.tolist()

    pas = pas.apply(lambda df: _pd_merge_gr(df,
                                            n_per_gene,
                                            how="left",
                                            on="gene_id",
                                            suffixes=["_pas","_counts"],
                                            to_merge_cols=per_gene_cols
                                            )
                    )

    #



if __name__ == '__main__':
    main()
