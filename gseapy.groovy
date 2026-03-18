container "quay.io/biocontainers/gseapy:1.1.11--py311h5e00ca1_1"

label "gseapy"
label "memory_medium"
label "process_medium"

when:((task.ext.when == null) || task.ext.when) && (de_file != null && de_file.size() > 0)

script:

run_gseapy = true //* @checkbox @description:"Run gene enrichment analysis using gseapy."

id_column = "" //* @input @description:"Name of the column containing gene IDs. Use 'index' if your table's first column is its index and it contains the gene IDs." @title:"Input file settings"
duplicate_id_strategy = "error" //* @dropdown @options:"error,keep_first,keep_last,remove_all" @description:"Strategy to handle duplicate gene IDs in the input gene set. 'error' will raise an error if duplicates are found. 'keep_first' will keep the first occurrence of each duplicate ID. 'keep_last' will keep the last occurrence of each duplicate ID. 'remove_all' will remove all occurrences of duplicate IDs."
pvalue_column = "" //* @input @description:"Name of the p-value column if the input gene set is to be filtered for a minimum p-value."
pvalue_cutoff = 0.05 //* @input @description:"Minimum p-value to filter the input gene set, only if pvalue_column is given."
lfc_column = "" //* @input @description:"Name of the fold change column if the input gene set is to be filtered for a minimum fold change."
lfc_cutoff = 0.0 //* @input @description:"Minimum log fold change to filter the input gene set, only if lfc_column is given."
rank_column = "stat" //* @input @description:"Name of the column to rank genes by for GSEA. Can be a statistic like t-statistic or log fold change. The column values should be numeric."
ascending = "no" //* @dropdown @options:"yes,no" @description:"Whether to rank genes in ascending order or not."
convert_to_uppercase = "yes" //* @dropdown @options:"yes,no" @description:"Whether to convert gene IDs to uppercase before enrichment analysis. This is useful to be able to use mouse genes with enrichr libraries that only support upper case even for mouse gene sets, as well as sources."
min_set_size = 15 //* @input @description:"Minimum gene set size for enrichment analysis."
max_set_size = 500 //* @input @description:"Maximum gene set size for enrichment analysis."
permutation_num = 1000 //* @input @description:"Number of permutations to perform for GSEA analysis. Higher numbers will give more accurate p-values but will take longer to run."

//* @style @condition:{run_gseapy=true,id_column,duplicate_id_strategy,pvalue_column,pvalue_cutoff,lfc_column,lfc_cutoff,organism,rank_column,ascending,msigdb_version,msig_dbs,enrichr_dbs,gmt_file,excel_file,convert_to_uppercase,min_set_size,max_set_size,permutation_num},{run_gprofiler2=false} @multicolumn:{id_column,duplicate_id_strategy,pvalue_column,pvalue_cutoff,lfc_column,lfc_cutoff,rank_column,ascending},{organism,msigdb_version,msig_dbs,enrichr_dbs,gmt_file,excel_file},{convert_to_uppercase,min_set_size,max_set_size,permutation_num}

id_column = task.ext.id_column ?: id_column
duplicate_id_strategy = task.ext.duplicate_id_strategy ?: duplicate_id_strategy
pvalue_column = task.ext.pvalue_column ?: pvalue_column
pvalue_cutoff = task.ext.pvalue_cutoff ?: pvalue_cutoff
lfc_column = task.ext.lfc_column ?: lfc_column
lfc_cutoff = task.ext.lfc_cutoff ?: lfc_cutoff
rank_column = task.ext.rank_column ?: rank_column
ascending = task.ext.ascending ?: ascending
convert_to_uppercase = task.ext.convert_to_uppercase ?: convert_to_uppercase
min_set_size = task.ext.min_set_size ?: min_set_size
max_set_size = task.ext.max_set_size ?: max_set_size
permutation_num = task.ext.permutation_num ?: permutation_num

query_name = de_file.simpleName

"""
#! /usr/bin/env python3
import os
import pandas as pd
import numpy as np
import gseapy as gp
import pickle
from gseapy import Msigdb as msig
import subprocess
from platform import python_version

de_file = "${de_file}"
library_file = "${library_file}"
id_column = "${id_column}"
duplicate_id_strategy = "${duplicate_id_strategy}"
pvalue_column = "${pvalue_column}"
pvalue_cutoff = ${pvalue_cutoff}
lfc_column = "${lfc_column}"
lfc_cutoff = ${lfc_cutoff}
query_name = "${query_name}"
rank_column = "${rank_column}"
ascending = "${ascending}" == "yes"
uppercase = "${convert_to_uppercase}" == "yes"
min_size = ${min_set_size}
max_size = ${max_set_size}
permutation_num = ${permutation_num}
threads = ${task.cpus}

# load library dictionary from file
with open(library_file, "rb") as f:
    library_dict = pickle.load(f)

# remove gene sets that are too small or too large
filtered_library_dict = {}
for lib_name, lib in library_dict.items():
    for gn, gene_sets in lib.items():
        fgs = {k:v for k,v in gene_sets.items()
               if len(v) >= min_size and len(v) <= max_size}

        # convert gene ids to uppercase if the option is set
        if uppercase:
            filtered_gene_sets = {k:list(map(str.upper, v)) for k,v in fgs.items()}
        else:
            filtered_gene_sets = fgs
        if filtered_gene_sets:
            try:
                filtered_library_dict[lib_name][gn] = filtered_gene_sets
            except KeyError:
                filtered_library_dict[lib_name] = {gn: filtered_gene_sets}
library_dict = filtered_library_dict

# read the DE results file
if id_column == "index":
    de = pd.read_csv(de_file, index_col=0)
    de.index.name = "gene"
    de.reset_index(inplace=True)
    id_column = "gene"
else:
    de = pd.read_csv(de_file)

# convert gene ids to uppercase if the option is set
if uppercase:
    de[id_column] = de[id_column].str.upper()

# check and handle duplicated gene ids
if de.duplicated(subset=id_column).any():
    dups = list(set(de.loc[de.duplicated(subset=id_column), id_column]))
    dups = "\\n".join(dups)
    if duplicate_id_strategy == "error":
        raise ValueError("Duplicated gene ids: {}".format(dups))
    else:
        keep_dict = {"remove_all": False,
                "keep_first": "first",
                "keep_last": "last"}
        drop_index = de.duplicated(subset=id_column,
                                   keep=keep_dict[duplicate_id_strategy])
        print("Duplicated gene IDs in the input: {}".format(dups))
        de = de.loc[~drop_index]

####### GSEA ANALYSIS WITH GSEAPY

# create ranked gene list
rnk = de.set_index(id_column).sort_values(rank_column, ascending=ascending)[
    rank_column]

# save the ranked gene list for potential future use
rnk.to_csv(query_name + "_ranked_gene_list.rnk", sep="\\t", header=False)

# run preranked GSEA for each gene set database
gsea_results = []
for lib_name, lib in library_dict.items():
    if len(lib) > 0:
        for gn, gene_sets in lib.items():
            if len(gene_sets) > 0:
                try:
                    pre_res = gp.prerank(rnk=rnk,
                                        gene_sets=gene_sets,
                                        min_size=min_size,
                                        max_size=max_size,
                                        threads=threads,
                                        permutation_num=permutation_num,
                                        outdir=None,
                                        verbose=True)

                    res= pre_res.res2d.copy()
                    res["Source"] = lib_name
                    res["Library"] = gn
                    gsea_results.append(res)
                except Exception as e:
                    print("Error running GSEA for library {}: {}".format(gn, str(e)))

if len(gsea_results) > 0:
    gsea_results = pd.concat(gsea_results, ignore_index=True)

    # add -log10 adjusted p-value and p-value to the results
    gsea_results["-logPadj"] = gsea_results["FDR q-val"].replace(
        0, gsea_results["FDR q-val"].replace(0, np.nan).min()).astype(float)
    gsea_results["-logPadj"] = -np.log10(gsea_results["-logPadj"])
    gsea_results["-logP"] = gsea_results["NOM p-val"].replace(
        0, gsea_results["NOM p-val"].replace(0, np.nan).min()).astype(float)
    gsea_results["-logP"] = -np.log10(gsea_results["-logP"])

    # add the query name to the results
    gsea_results["Query"] = query_name

    # add positive or negative enrichment info
    gsea_results["Enrichment"] = gsea_results["NES"].apply(
        lambda a: "Positive" if a > 0 else "Negative")

    # process tag% column to get set size and other useful values
    gsea_results[["Overlap Size", "Term Size"]] = gsea_results["Tag %"].apply(
        lambda a: pd.Series(map(int, a.split("/"))))
    gsea_results["Overlap %"] = round(
        gsea_results["Overlap Size"] / gsea_results["Term Size"] * 100, 1)

    # convert gene % values to numbers
    gsea_results["Gene %"] = gsea_results["Gene %"].apply(
        lambda a: round(float(a.replace("%", "")), 1))
    # save gsea results
    gsea_results.to_csv(query_name + "_gsea_results.csv", index=False)

else:
    with open(query_name + "_gsea_results.csv", "w") as f:
        pass

#### ORA ANALYSIS WITH GSEAPY

# prepare DE gene list for ORA
# filter significant genes based on p-value and log fold change cutoffs
if pvalue_column not in ("", "NA", "na"):
    sig_mask = de[pvalue_column] <= pvalue_cutoff
    sig_df = de.loc[sig_mask].copy()
else:
    sig_df = de.copy()

enrich_results_filename = query_name + "_enrich_results.csv"

# create empty results file if no significant genes are found after filtering
if sig_df.empty:
    with open(enrich_results_filename, "w") as f:
        pass
else:
    if lfc_column not in ("", "NA", "na"):
        up_mask = de[lfc_column] >= lfc_cutoff
        down_mask = de[lfc_column] <= lfc_cutoff
        up_genes = list(sig_df.loc[up_mask, id_column])
        down_genes = list(sig_df.loc[down_mask, id_column])
        lfc_provided = True
    else:
        lfc_provided = False
        sig_genes = list(sig_df[id_column])

    # create query dictionaries
    if not lfc_provided:
        q_dict = {query_name: sig_genes}
    else:
        q_dict = {}
        if up_genes:
            q_dict[query_name + "_UP"] = up_genes
        if down_genes:
            q_dict[query_name + "_DN"] = down_genes

    # use all genes as background
    background_genes = list(de[id_column])

    # run ORA for each gene set database
    enrich_results = []
    for lib_name, lib in library_dict.items():
        for gn, gene_sets in lib.items():
            if len(gene_sets) > 0:
                for q, gene_list in q_dict.items():
                    enrich_res = gp.enrich( gene_list=gene_list,
                            gene_sets=gene_sets,
                            background=background_genes,
                            no_plot=True,
                            outdir=None,
                        verbose=False)
                    res = enrich_res.results.copy()

                    if len(res) > 0:
                        res= enrich_res.results.copy()
                        res["Source"] = lib_name
                        res["Library"] = gn
                        res["Query"] = q
                        enrich_results.append(res)

    if len(enrich_results) > 0:
        enrich_results = pd.concat(enrich_results, ignore_index=True)

        # add -log10 adjusted p-value to the results
        enrich_results["-logPadj"] = -np.log10(enrich_results["Adjusted P-value"])
        enrich_results["-logP"] = -np.log10(enrich_results["P-value"])

        # add a nicer fomatted version of up or down regulated query names
        if lfc_provided:
            enrich_results["Expression"] = enrich_results["Query"].apply(
                lambda a: "Up-regulated" if a.split("_")[-1] == "UP"
                    else "Down-regulated")

        # process overlap column to get set size and other useful values
            enrich_results[["Overlap Size", "Term Size"]] = enrich_results["Overlap"].apply(
                lambda a: pd.Series(map(int, a.split("/"))))
            enrich_results["Overlap %"] = round(
                enrich_results["Overlap Size"] / enrich_results["Term Size"] * 100, 1)

        # save results
        enrich_results.to_csv(enrich_results_filename, index=False)

    else:
        with open(enrich_results_filename, "w") as f:
            pass

# versions

versions = {"python": python_version(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "gseapy": gp.__version__,
            "msigdb": msigdb_version}

with open("versions.yml", "w") as outfile:
    outfile.write("$task.process" + ":\\n")
    for v in versions:
        outfile.write("\\t" + v + ": " + versions[v] + "\\n")

subprocess.call(["cp", ".command.sh", query_name + ".${task.process}.command.sh"])
subprocess.call(["cp", ".command.err", query_name + ".gseapy.log"])

"""
