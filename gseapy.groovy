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
organism = "" //* @dropdown @options:"mouse,human,yeast,worm,fly,fish" @description:"Organism name. Limited support for yeast, worm, fly and fish."
msigdb_version = "default" //* @input @description:"Version of the MSigDB gene sets to use. Default is 2026.1.Mm for mouse and 2026.1.Hs for human."
msig_dbs = "default" //* @input @description:"Comma-separated list of MSigDB gene set databases to use. Default is a set of selected databases for mouse and human. See pipeline documentation for options."
enrichr_dbs = "default" //* @input @description:"Comma-separated list of Enrichr gene set libraries to use. Default is a set of selected libraries. See https://maayanlab.cloud/Enrichr/#libraries for options."
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
organism = task.ext.organism ?: organism
rank_column = task.ext.rank_column ?: rank_column
ascending = task.ext.ascending ?: ascending
msigdb_version = task.ext.msigdb_version ?: msigdb_version
msig_dbs = task.ext.msig_dbs ?: msig_dbs
enrichr_dbs = task.ext.enrichr_dbs ?: enrichr_dbs
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
from gseapy import Msigdb as msig
import subprocess
from platform import python_version

de_file = "${de_file}"
id_column = "${id_column}"
duplicate_id_strategy = "${duplicate_id_strategy}"
pvalue_column = "${pvalue_column}"
pvalue_cutoff = ${pvalue_cutoff}
lfc_column = "${lfc_column}"
lfc_cutoff = ${lfc_cutoff}
query_name = "${query_name}"
organism = "${organism}"
rank_column = "${rank_column}"
ascending = "${ascending}" == "yes"
msigdb_version = "${msigdb_version}"
msig_dbs = "${msig_dbs}"
enrichr_dbs = "${enrichr_dbs}"
gmt_file = "${gmt_file}"
excel_file = "${excel_file}"
uppercase = "${convert_to_uppercase}" == "yes"
min_size = ${min_set_size}
max_size = ${max_set_size}
permutation_num = ${permutation_num}
threads = ${task.cpus}

# create a dictionary to hold all gene set libraries
library_dict = {}

# get gene sets from input sources
# 1. msigdb gene sets

# set default MSigDB gene sets if not provided by user
human_msig_dbs = (
    'h.all, c2.cgp, c2.cp.biocarta, c2.cp.kegg_medicus, c2.cp.pid,'
    'c2.cp.reactome, c2.cp.wikipathways, c3.mir.mirdb, c3.tft.gtrd, c5.go.bp, '
    'c5.go.cc, c5.go.mf, c7.immunesigdb, c9.all'
)
mouse_msig_dbs = (
    'mh.all, m2.cgp, m2.cp.biocarta, m2.cp.reactome, m2.cp.wikipathways, m3.gtrd,'
    'm3.mirdb, m5.go.bp, m5.go.cc, m5.go.mf, m7.all'
)

# set test MSigDB gene sets for testing purposes
human_test_msig_dbs = "c2.cp.reactome, c5.go.bp"
mouse_test_msig_dbs = "m2.cp.reactome, m5.go.bp"

if msig_dbs == "default":
    if organism == "mouse":
        msig_dbs = mouse_msig_dbs
    elif organism == "human":
        msig_dbs = human_msig_dbs
elif msig_dbs == "test":
    if organism == "mouse":
        msig_dbs = mouse_test_msig_dbs
    elif organism == "human":
        msig_dbs = human_test_msig_dbs

# convert the comma-separated string of gene set databases into a list
msig_dbs = list(set([db.strip() for db in msig_dbs.strip().split(",")]))

# set default MSigDB version if not provided by user
if msigdb_version == "default":
    if organism == "mouse":
        msigdb_version = "2026.1.Mm"
    elif organism == "human":
        msigdb_version = "2026.1.Hs"

# get gene sets from MSigDB
msig_gs_dict = {}

for db_name in msig_dbs:
    db = msig.get_gmt(category=db_name, dbver=msigdb_version)
    if db is not None:
        if uppercase:
            db = {k:list(map(str.upper, v)) for k,v in db.items()}
        msig_gs_dict[db_name] = db
    else:
        print("MSig gene set {} returned no genes for database version {}".format(
            db_name, msigdb_version))
if msig_gs_dict:
    library_dict["msigdb"] = msig_gs_dict

# 2. Enrichr gene sets

# set default Enrichr gene set libraries if not provided by user
if enrichr_dbs == "default":
    enrichr_dbs = (
        "BioCarta_2016, BioPlanet_2019, CORUM, Elsevier_Pathway_Collection, "
        "KEGG_2026, ARCHS4_TFs_Coexp, ChEA_2022, "
        "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X,"
        "DepMap_CRISPR_GeneDependency_CellLines_2023, Panther_2016, "
        "HMDB_Metabolites, JASPAR_PWM_Human_2025, TRANSFAC_and_JASPAR_PWMs,"
        "TargetScan_microRNA_2017, ARCHS4_Kinases_Coexp, "
        "Kinase_Perturbations_from_GEO_down, Kinase_Perturbations_from_GEO_up, "
        "L1000_Kinase_and_GPCR_Perturbations_down,"
        "L1000_Kinase_and_GPCR_Perturbations_up"
    )
elif enrichr_dbs == "test":
    enrichr_dbs = "ChEA_2022, KEGG_2026"

enrichr_dbs = list(set([db.strip() for db in enrichr_dbs.strip().split(",")]))

# get gene sets from Enrichr
enrichr_gs_dict = {}

for db_name in enrichr_dbs:
    db = gp.get_library(name=db_name, organism=organism)
    if db is not None:
        enrichr_gs_dict[db_name] = db
    else:
        print("Enrichr library {} returned no genes for organism {}".format(
            db_name, organism))
if enrichr_gs_dict:
    library_dict["enrichr"] = enrichr_gs_dict

# 3. gene sets from user-provided excel file
try:
    gs_excel = pd.read_excel(excel_file, header=None, sheet_name=None)
except (ValueError, NameError):
    gs_excel = {}
excel_gene_sets = {}
for k,v in gs_excel.items():
    s_name = "".join([c if c.isalnum() else "_" for c in k.strip()])
    s_list = list(v.squeeze().str.strip().unique())
    if uppercase:
        s_list = list(map(str.upper, s_list))
    excel_gene_sets[s_name] = s_list
if excel_gene_sets:
    library_dict["excel"] = {"excel": excel_gene_sets}

# 4. gene sets from user-provided GMT file
try:
    gmt_gene_sets = gp.read_gmt(gmt_file)
    if uppercase:
        gmt_gene_sets = {k:list(map(str.upper, v)) for k,v in gmt_gene_sets.items()}
except (ValueError, NameError):
    gmt_gene_sets = {}
if gmt_gene_sets:
    library_dict["gmt"] = {"gmt": gmt_gene_sets}

# remove gene sets that are too small or too large
filtered_library_dict = {}
for lib_name, lib in library_dict.items():
    for gn, gene_sets in lib.items():
        filtered_gene_sets = {k:v for k,v in gene_sets.items()
                                if len(v) >= min_size and len(v) <= max_size}
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
