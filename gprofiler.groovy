container "quay.io/umassbio/requests:2.31.0--0"

label "requests"

when: ( task.ext.when == null ) || task.ext.when

script:
run_gprofiler = true //* @checkbox @description:"Run gene enrichment analysis using gprofiler."

id_column = "" //* @input @description:"Name of the column containing gene IDs. Use 'index' if your table's first column is its index and it contains the gene IDs." @title:"Input file settings"
pvalue_column = "" //* @input @description:"Name of the p-value column if the input gene set is to be filtered for a minimum p-value."
pvalue_cutoff = 0.05 //* @input @description:"Minimum p-value to filter the input gene set, only if pvalue_column is given."
lfc_column = "" //* @input @description:"Name of the fold change column if the input gene set is to be filtered for a minimum fold change."
lfc_cutoff = 0.0 //* @input @description:"Minimum log fold change to filter the input gene set, only if lfc_column is given."
universe = "no" //* @dropdown @options:"yes,no" @description:"By default, the background for enrichment uses all the genes in the input dataset. Set this to yes if you prefer using all (annotated or known) genes of the species even if the genes are not included in your data set." @title:"Statistics settings"
stats_domain = "annotated" //* @dropdown @options:"annotated,known" @description:"Use all genes in the background gene set (known) or only those that have annotations (annotated)."
significant_only = "yes" //* @dropdown @options:"yes,no" @description:"Report only significan results."
padj_method = "g_SCS" //* @dropdown @options:"g_SCS","bonferroni","fdr","analytical" @description:"Multiple testing correction method to use."
padj_cutoff = 0.05 //* @input @description:"Minimum adjusted p-value for output terms."
exclude_iea = "no" //* @dropdown @options:"yes,no" @description:"Exclude electronically annotated terms."
duplicate_id_strategy = "error" //* @dropdown @options:"error,keep_first,keep_last,remove_all" @description:"Strategy to handle duplicate gene IDs in the input gene set. 'error' will raise an error if duplicates are found. 'keep_first' will keep the first occurrence of each duplicate ID. 'keep_last' will keep the last occurrence of each duplicate ID. 'remove_all' will remove all occurrences of duplicate IDs."
organism = "" //* @input @description:"Organism name. Use the g:Profiler short name for the organism. For example, use 'hsapiens' for human and 'mmusculus' for mouse. See https://biit.cs.ut.ee/gprofiler/page/species for the full list of supported organisms."
//* @style @condition:{run_gprofiler=true,id_column,duplicate_id_strategy,pvalue_column,pvalue_cutoff,stats_filter,lfc_column,lfc_cutoff,organism,universe,stats_domain,significant_only,padj_method,padj_cutoff,exclude_iea},{run_gprofiler2=false} @multicolumn:{id_column,duplicate_id_strategy,pvalue_column,pvalue_cutoff,lfc_column,lfc_cutoff},{organism,universe,stats_domain,exclude_iea,significant_only,padj_method,padj_cutoff}

id_column = task.ext.id_column ?: id_column
pvalue_column = task.ext.pvalue_column ?: pvalue_column
lfc_column = task.ext.lfc_column ?: lfc_column
lfc_cutoff = task.ext.lfc_cutoff ?: lfc_cutoff
pvalue_cutoff = task.ext.pvalue_cutoff ?: pvalue_cutoff
universe = task.ext.universe ?: universe
stats_domain = task.ext.stats_domain ?: stats_domain
significant_only = task.ext.significant_only ?: significant_only
padj_method = task.ext.padj_method ?: padj_method
padj_cutoff = task.ext.padj_cutoff ?: padj_cutoff
exclude_iea = task.ext.exclude_iea ?: exclude_iea
duplicate_id_strategy = task.ext.duplicate_id_strategy ?: duplicate_id_strategy
organism = task.ext.organism ?: organism
organism = organism == "" ? params.gprofiler_name : organism
query_name = de_file.simpleName

"""
#! /usr/bin/env python3

import os
import pandas as pd
import numpy as np
import requests
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
stats_domain = "${stats_domain}"
universe = "${universe}"
significant_only = "${significant_only}"
padj_method = "${padj_method}"
padj_cutoff = ${padj_cutoff}
exclude_iea = "${exclude_iea}"

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

# filter significant genes based on p-value and log fold change cutoffs   
if pvalue_column not in ("", "NA", "na"):
    sig_mask = de[pvalue_column] <= pvalue_cutoff
    sig_df = de.loc[sig_mask].copy()
else:
    sig_df = de.copy()
    
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
    q_dict = {query_name + "_UP": up_genes,
              query_name + "_DN": down_genes}

# get all genes as background if universe is not provided
if universe == "no":
    stats_domain = "custom_" + stats_domain
    background_genes = list(de[id_column])
else:
    background_genes = None
    
# create the payload for g:Profiler API
query_dict = {'organism': organism,
              'query': q_dict,
              "background": background_genes,
              'sources' :[],
              'domain_scope': stats_domain,
              "highlight": significant_only == "yes",
              "all_results": significant_only != "yes",
              "user_threshold": padj_cutoff,
              "significance_threshold_method": padj_method,
              "no_iea": exclude_iea == "yes",
              "no_evidences": False}

# send the request to g:Profiler API
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
    json=query_dict,
    headers={'User-Agent':'FullPythonRequest'}
)    

# check the response status
if not r.ok:
    msg= (("API request failed with status code: {}. "
           "Message from the gprofiler server was: {}")).format(r.status_code, r.text)
    raise Exception(msg)

# get the results as a dataframe
res = r.json()
results = pd.DataFrame.from_records(res["result"])

# get the metadata for the query
meta = res["meta"]

# create ensemble gene id to gene name mapping for queries
for qry, qry_dict in meta["genes_metadata"]["query"].items():
    ens = qry_dict["ensgs"]
    qry_dict["ensgs_array"] = np.array(ens)
    qry_dict["ensgs_to_gene"] = {}
    gene_to_ens = qry_dict["mapping"]
    for gene_name, ensembl_list in gene_to_ens.items():
        for ens_id in ensembl_list:
            try:
               qry_dict["ensgs_to_gene"][ens_id].update([gene_name])
            except KeyError:
                qry_dict["ensgs_to_gene"][ens_id] = set([gene_name])

# define a function to convert ensembl gene ids to gene names in the results
def get_intersection_genes(row):
    q = row["query"]
    # get ensembl IDs of all genes in query
    ens_array = meta["genes_metadata"]["query"][q]["ensgs_array"]
    # get indexes of genes that are intersecting the term
    ind = np.array([(len(i) > 0) for i in row["intersections"]])
    # get ensembl IDs of genes in the intersection
    ens = ens_array[ind]
    # get ensembl to gene mapping for query
    mapping_dict = meta["genes_metadata"]["query"][q]["ensgs_to_gene"]
    # create a set to add gene ids matching the ensembl IDs in the intersection
    genes = set()
    for e in ens:
        genes.update(mapping_dict[e])
    return genes

# apply the function to get the gene names for the intersecting genes
results["genes"] = results.apply(get_intersection_genes, axis=1)

# calcutale fold enrichment for each term
results["Fold Enrichment"] = (
    (results["intersection_size"] / results["query_size"])
    /(results["term_size"] / results["effective_domain_size"]))

# add database source names to the results
db_map = {"GO:CC": "GO Cellular Component",
          "GO:BP": "GO Biological Process",
          "GO:MF": "GO Molecular Function",
          "HP": "Human Protein Atlas",
          "WP": "WikiPathways",
          "KEGG": "KEGG",
          "TF": "TRANSFAC",
          "MIRNA": "miRTarBase",
          "CORUM": "CORUM",
          "REAC": "Reactome"}

results["database"] = results["source"].replace(db_map)

# add -log10 adjusted p-value to the results
results["-logPadj"] = -np.log10(results["p_value"])

# add a nicer fomatted version of up or down regulated query names
if lfc_provided:
    results["Expression"] = results["query"].apply(
        lambda a: "Up-regulated" if a.split("_")[-1] == "UP"
              else "Down-regulated")

# save the results to a csv file
results.drop(columns={"intersections"}, inplace=True)
results.to_csv(query_name + ".gprofiler.csv", index=False)

# versions

versions = {"python": python_version(),
			"pandas": pd.__version__,
            "numpy": np.__version__,
            "requests": requests.__version__,
			"gprofiler_api": meta["version"],
			"gprofiler_timestamp": meta["timestamp"]
}
with open("versions.yml", "w") as outfile:
    outfile.write("$task.process" + ":\\n")
    for v in versions:
    	outfile.write("\\t" + v + ": " + versions[v] + "\\n")

subprocess.call(["cp", ".command.sh", query_name + ".${task.process}.command.sh"])
"""
