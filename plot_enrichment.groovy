
label "gseapy"

when: ((task.ext.when == null) || task.ext.when) && (enrichment_file != null && enrichment_file.size() > 0)

script:

x_column = "" //* @input @description:"Column name for the x-axis values." @title:"Input file settings"
y_column = "" //* @input @description:"Column name for the y-axis values."
size_column = "" //* @input @description:"Column name for the size of points in a dotplot."
min_point_size = 10 //* @input @description:"Minimum point size for dotplot."
max_point_size = 200 //* @input @description:"Maximum point size for dotplot."
split_columns = "" //* @input @description:"Comma-separated list of column names to split the input data by for plotting. Separate plots will be generated for each unique combination of values in these columns."
plot_title_column = "" //* @input @description:"Column name to use for plot titles when split_columns is specified. If not provided, last column in split_columns will be used for plot titles."
group_column = "" //* @input @description:"Column name for grouping variable to color bars or points by."
group_order = "" //* @input @description:"Comma-separated list of group names in the order they should appear in the plot. Only applicable if group_column is specified."
group_colors = "" //* @input @description:"Comma-separated list of colors to use for each group. Only applicable if group_column is specified. If not provided, default colors will be used."
single_color = "" //* @input @description:"Color to use for bars or points when group_column is not specified."
sort_ascending = "no" //* @dropdown @options:"yes,no" @description:"Whether to sort the data in ascending order based on x_column for selecting top N results."
top_n = 10 //* @input @description:"Number of top results to show in the plot. If group_column is specified, top N results will be selected for each group. Use -1 to show all results (maximum 50)."")"
pvalue_column = "" //* @input @description:"Name of the p-value column if the input is to be filtered for a minimum p-value."
pvalue_cutoff = 0.05 //* @input @description:"Minimum p-value to filter the input, only if pvalue_column is given."

//* @style @multicolumn:{x_column, y_column, size_column, group_column, group_order, group_colors, single_color},{sort_ascending, top_n, pvalue_column, pvalue_cutoff}

x_column = task.ext.x_column ?: x_column
y_column = task.ext.y_column ?: y_column
size_column = task.ext.size_column ?: size_column
min_point_size = task.ext.min_point_size ?: min_point_size
max_point_size = task.ext.max_point_size ?: max_point_size
split_columns = task.ext.split_columns ?: split_columns
plot_title_column = task.ext.plot_title_column ?: plot_title_column
group_column = task.ext.group_column ?: group_column
group_order = task.ext.group_order ?: group_order
group_colors = task.ext.group_colors ?: group_colors
single_color = task.ext.single_color ?: single_color
sort_ascending = task.ext.sort_ascending ?: sort_ascending
top_n = task.ext.top_n ?: top_n
pvalue_column = task.ext.pvalue_column ?: pvalue_column
pvalue_cutoff = task.ext.pvalue_cutoff ?: pvalue_cutoff

"""
!# /usr/bin/env python3
import matplotlib.pyplot as plt # type: ignore
import seaborn as sns
import pandas as pd # type: ignore
import numpy as np # type: ignore
import subprocess
from platform import python_version

# convert nextflow parameters to python variables
x_column = "!{x_column}"
y_column = "!{y_column}"
size_column = "!{size_column}" if "!{size_column}" else None
min_point_size = !{min_point_size}
max_point_size = !{max_point_size}
dot_sizes = (min_point_size, max_point_size)
split_columns = [c.strip() for c in "!{split_columns}".split(",") if c.strip()]
plot_title_column = "!{plot_title_column}"
group_column = "!{group_column}" if "!{group_column}" else None
group_order = [g.strip() for g in "!{group_order}".split(",") if g.strip()]
if not group_order:
	group_order = None
group_colors = [c.strip() for c in "!{group_colors}".split(",") if c.strip()]
if not group_colors:
	group_colors = None
single_color = "!{single_color}" if "!{single_color}" else None
sort_ascending = "!{sort_ascending}".lower() == "yes"
top_n = int(!{top_n})
if top_n == -1:
	top_n = None
pvalue_column = "!{pvalue_column}" if "!{pvalue_column}" else None
pvalue_cutoff = !{pvalue_cutoff}
query_name = "!{query_name}"
enrichment_file = "!{enrichment_file}"

# define a function to create enrichment plots
def enrichment_plot(input_df, x_col, y_col, plot_type,
                    p_col=None, p_cut=None, size_col=None, dot_sizes=(10, 200),
                    group_col=None, group_order=None, group_colors=None,
                    color=None, sort_ascending=False, top_n=10,
                    plot_title=None):
    """
    Plot enrichment results as 'barplot' or 'dotplot'.
    """

    # work on a copy of the input data to avoid modifying the caller's DataFrame
	df = input_df.copy()

    # filter significant results if p_cut is specified
    if p_col is not None and p_cut is not None:
        if p_col not in df.columns:
            raise KeyError(f"p_col '{p_col}' not found in enrichment_results")
        sig_results = df[df[p_col] < p_cut].copy()
    else:
        sig_results = df

    # ensure x_col and y_col exist in the data to plot
    for col in (x_col, y_col):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in the data to plot")

    # set group color palette to None (may be created later)
    palette = None

    # select top N results for each group (no extra index levels)
    if top_n is None:
        top_n = len(sig_results)

    if group_col is not None:
        if group_col not in sig_results.columns:
            raise KeyError(f"group_col '{group_col}' not found in enrichment_results")

        top_data = (
            sig_results
            .groupby(group_col)
            .apply(lambda a: a.sort_values(x_col, ascending=sort_ascending).head(top_n))
            .reset_index()
        )

        # add a newline to the beginning of the size col to provide
        # some separation between the hue and size legends
        if size_col is not None and size_col in top_data.columns:
            new_size_col = "\n" + size_col
            top_data = top_data.rename(columns={size_col: new_size_col})
            size_col = new_size_col

        # sort by group and x_col if group order is specified
        if group_order is not None:
            order_map = {g: i for i, g in enumerate(group_order)}
            top_data["group_order"] = top_data[group_col].map(order_map)
            top_data = top_data.sort_values(["group_order", x_col],
                                            ascending=[True, sort_ascending]
                                            ).drop(columns=["group_order"])
            # create color palette if group colors are specified
            if group_colors is not None:
                if len(group_colors) != len(group_order):
                    raise ValueError("Length of group_colors must match length of group_order")
                palette = dict(zip(group_order, group_colors))

    else:
        top_data = sig_results.sort_values(
            x_col, ascending=sort_ascending).head(top_n).reset_index(drop=True)
        # reset group related parameters in case they are provided but not used
        group_order = None
        group_colors = None

    if plot_type == "barplot":
        fig, ax = plt.subplots()

        # build kwargs to avoid passing color when hue/group is used
        bar_kwargs = dict(data=top_data, y=y_col, x=x_col, orient="h",
                          hue=group_col, hue_order=group_order, palette=palette,
                          dodge=False, ax=ax)
        if color is not None and group_col is None:
            bar_kwargs["color"] = color

        p = sns.barplot(**bar_kwargs)
        fig.set_size_inches(5, max(15, top_data.shape[0] / 5))
        ax.set_ylabel("")
        ax.set_title(plot_title)

        try:
            # move legend using the axes object
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        except Exception:
            pass

        return fig, ax

    elif plot_type == "dotplot":
        with sns.axes_style("whitegrid"):
            fig, ax = plt.subplots()

            scatter_kwargs = dict(data=top_data, x=x_col, y=y_col,
                                  size=size_col, sizes=(10, 200),
                                  ax=ax, hue=group_col, legend="brief",
                                  palette=palette, hue_order=group_order)
            if color is not None and group_col is None:
                scatter_kwargs["color"] = color

            sns.scatterplot(**scatter_kwargs)
            ax.grid(visible=False, axis="x")
            ax.set_ylabel("")
            ax.set_title(plot_title)

            try:
                sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
            except Exception:
                pass

            return fig, ax

    else:
        raise ValueError(f"Unknown plot type: {plot_type}. Supported plot types are: ['barplot', 'dotplot'].")


# load enrichment results from file
enrichment_results = pd.read_csv("${enrichment_file}")

# split data by specified columns and generate a plot for each subset
enrichment_results["temp_col"] = " "
if split_columns:
	plot_title_column = split_columns[-1]
else:
	split_columns = ["temp_col"]
	if plot_title_column:
		if len(enrichment_results[plot_title_column].unique()) > 1:
			raise ValueError(
			f"Error: plot_title_column '{plot_title_column}' has multiple unique "
			"values but split_columns is not specified. When split_columns is not "
			"specified, title column must have only one unique value to be used for"
			" plot titles.")
		else:
			plot_title_column = "temp_col"

grouped = enrichment_results.groupby(split_columns, as_index=False)
for group_vals, group_data in grouped:
	plot_title = query_name + "\n" + group_data[plot_title_column].iloc[0]
	split_vals = '_'.join([v.strip() for v in group_vals])
	file_prefix = f"{query_name}_{split_vals}"
	fig, ax = enrichment_plot(group_data, x_column, y_column, "dotplot",
							size_col=size_column, dot_sizes=dot_sizes,
							group_col=group_column,
							group_order=group_order, group_colors=group_colors,
							color=single_color, sort_ascending=sort_ascending,
							top_n=top_n, plot_title=plot_title)
	fig.savefig(f"{file_prefix}_enrichment_dotplot.pdf", bbox_inches="tight")

	fig, ax = enrichment_plot(group_data, x_column, y_column, "barplot",
							group_col=group_column, group_order=group_order,
							group_colors=group_colors, color=single_color,
							sort_ascending=sort_ascending, top_n=top_n,
							plot_title=plot_title)
	fig.savefig(f"{file_prefix}_enrichment_barplot.pdf", bbox_inches="tight")


# versions

versions = {"python": python_version(),
			"pandas": pd.__version__,
            "numpy": np.__version__,
			"seaborn": sns.__version__,
			"matplotlib": plt.matplotlib.__version__}

with open("versions.yml", "w") as outfile:
    outfile.write("$task.process" + ":\\n")
    for v in versions:
    	outfile.write("\\t" + v + ": " + versions[v] + "\\n")

subprocess.call(["cp", ".command.sh", query_name + ".${task.process}.command.sh"])
"""
