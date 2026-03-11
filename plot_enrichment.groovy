
label "seaborn"

when: ((task.ext.when == null) || task.ext.when) && (enrichment_file != null && enrichment_file.size() > 0)

script:

plot_type = "dotplot" //* @dropdown @options:"dotplot,barplot" @description:"Type of plot to generate for enrichment results."
x_column = "" //* @input @description:"Column name for the x-axis values." @title:"Plot settings"
y_column = "" //* @input @description:"Column name for the y-axis values."

split_columns = "" //* @input @description:"Comma-separated list of column names to split the input data by for plotting. Separate plots will be generated for each unique combination of values in these columns." @title:"Split/group settings"
plot_title_column = "" //* @input @description:"Column name to use for plot titles when split_columns is not specified. This column must have a single unique value."
group_column = "" //* @input @description:"Column name for grouping within the same plot. For example, up-regulated vs down-regulated."
group_order = "" //* @input @description:"Comma-separated list of group names in the order they should appear in the plot. Only applicable if group_column is specified."
group_colors = "" //* @input @description:"Comma-separated list of colors to use for each group. Only applicable if group_column is specified. If not provided, default colors will be used."

hue_column = "" //* @input @description:"Column name for variable to color points by in dotplot. Used only if group_column is not specified."  @title:"Color settings"
hue_order = "" //* @input @description:"Comma-separated list of hue variable values in the order they should appear in the plot. Only applicable if hue_column is specified."
hue_values = "" //* @input @description:"Comma-separated list of colors to use for each hue variable value."
palette = "" //* @input @description:"Name of matplotlib/seaborn color palette to use. Only applicable if group_column or hue_column is specified and group_colors or hue_values are not provided. Single color sequential palettes (Greys,Reds, Blues, Oranges, Greens, Purples) are recommended for quantitative hue columns, such as p-values or overlap %. See https://seaborn.pydata.org/tutorial/color_palettes.html for more details on available color palettes."
single_color = "" //* @input @description:"Color to use for bars or points when group_column is not specified."

sort_ascending = "no" //* @dropdown @options:"yes,no" @description:"Whether to sort the data in ascending order based on x_column for selecting top N results." @title:"Misc settings"
top_n = 10 //* @input @description:"Number of top results to show in the plot. If group_column is specified, top N results will be selected for each group. Use -1 to show all results (maximum 50)."")"
fig_size = "" //* @input @description:"Width and height of the plot in inches, separated by a comma (e.g. 5,7). If not provided, size will be calculated based on the number of data points."
pvalue_column = "" //* @input @description:"Name of the p-value column if the input is to be filtered for a minimum p-value."
pvalue_cutoff = 0.05 //* @input @description:"Minimum p-value to filter the input, only if pvalue_column is given."

size_column = "" //* @input @description:"Column name for the size of points in a dotplot." @title:"Dot plot settings"
min_point_size = 5 //* @input @description:"Minimum point size for dotplot."
max_point_size = 15 //* @input @description:"Maximum point size for dotplot."
//* @style @multicolumn: {plot_type, x_column, y_column}, {split_columns, plot_title_column, group_column, group_order, group_colors}, {hue_column, hue_order, hue_values, palette, single_color}, {sort_ascending, top_n, fig_size, pvalue_column, pvalue_cutoff}, {size_column, min_point_size, max_point_size}

plot_type = task.ext.plot_type ?: plot_type
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
hue_column = task.ext.hue_column ?: hue_column
hue_order = task.ext.hue_order ?: hue_order
hue_values = task.ext.hue_values ?: hue_values
palette = task.ext.palette ?: palette
single_color = task.ext.single_color ?: single_color
sort_ascending = task.ext.sort_ascending ?: sort_ascending
top_n = task.ext.top_n ?: top_n
fig_size = task.ext.fig_size ?: fig_size
pvalue_column = task.ext.pvalue_column ?: pvalue_column
pvalue_cutoff = task.ext.pvalue_cutoff ?: pvalue_cutoff

query_name = task.ext.query_name ?: query_name

"""
!# /usr/bin/env python3
import matplotlib.pyplot as plt # type: ignore
import seaborn as sns
from seaborn import axes_style
import seaborn.objects as so
import pandas as pd # type: ignore
from pandas.api.types import is_numeric_dtype
import numpy as np # type: ignore
from cycler import cycler
import subprocess
from platform import python_version

# convert nextflow parameters to python variables
plot_type = "!{plot_type}"
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

hue_column = "!{hue_column}" if "!{hue_column}" else None
hue_order = [h.strip() for h in "!{hue_order}".split(",") if h.strip()]
if not hue_order:
    hue_order = None
hue_values = [c.strip() for c in "!{hue_values}".split(",") if c.strip()]
if not hue_values:
    hue_values = None

palette = "!{palette}" if "!{palette}" else None
single_color = "!{single_color}" if "!{single_color}" else None
sort_ascending = "!{sort_ascending}".lower() == "yes"
top_n = int(!{top_n})
if top_n == -1:
    top_n = None
pvalue_column = "!{pvalue_column}" if "!{pvalue_column}" else None
pvalue_cutoff = !{pvalue_cutoff}
fig_size = "!{fig_size}" if "!{fig_size}" else None
if fig_size:
    try:
        fig_size = [float(c.strip()) for c in fig_size.split(",") if c.strip()]
        if len(fig_size) != 2:
            raise ValueError("fig_size must have two comma-separated values for width and height in inches.")
    except Exception:
        raise ValueError(f"fig_size: {fig_size} must have two comma-separated values for width and height in inches.")
query_name = "!{query_name}"
query_name = "" if query_name.lower() in ["none", "na", "null", ""] else query_name

enrichment_file = "!{enrichment_file}"

# define a function to create enrichment plots
def enrichment_plot(input_df, x_column, y_column, plot_type,
                    pvalue_column=None, pvalue_cutoff=None, size_column=None, dot_sizes=(5, 20),
                    group_column=None, group_order=None, group_colors=None,
                    hue_column=None, hue_order=None, hue_values=None,
                    palette=None, single_color=None, sort_ascending=False, top_n=10,
                    plot_title=None, fig_size=None):

    # work on a copy of the input data to avoid modifying the caller's DataFrame
    df = input_df.copy()

    # filter significant results if pvalue_cutoff is specified
    if pvalue_column is not None and pvalue_cutoff is not None:
        if pvalue_column not in df.columns:
            raise KeyError(f"pvalue_column '{pvalue_column}' not found in enrichment_results")
        sig_results = df[df[pvalue_column] < pvalue_cutoff].copy()
    else:
        sig_results = df

    # ensure x_column and y_column exist in the data to plot
    for col in (x_column, y_column):
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in the data to plot")

    # select top N results for each group (no extra index levels)
    if top_n is None:
        top_n = len(sig_results)

    nominal_color = True
    if group_column is not None:
        if group_column not in sig_results.columns:
            raise KeyError(f"group_column '{group_column}' not found in enrichment_results")

        top_data = (
            sig_results
            .groupby(group_column)
            .apply(lambda a: a.sort_values(x_column, ascending=sort_ascending).head(top_n))
            .reset_index()
        )

        # sort by group and x_column if group order is specified
        if group_order is not None:
            order_map = {g: i for i, g in enumerate(group_order)}
            top_data["group_order"] = top_data[group_column].map(order_map)
            top_data = top_data.sort_values(["group_order", x_column],
                                            ascending=[True, sort_ascending]
                                            ).drop(columns=["group_order"])

            # create color palette if group colors are specified
            if group_colors is not None:
                if len(group_colors) != len(group_order):
                    raise ValueError("Length of group_colors must match length of group_order")
                palette = dict(zip(group_order, group_colors))

        hue_column = group_column
        hue_order = group_order

    else:
        top_data = sig_results.sort_values(
            x_column, ascending=sort_ascending).head(top_n).reset_index(drop=True)
        if hue_column is not None:
            if hue_column not in sig_results.columns:
                raise KeyError(f"group_column '{group_column}' not found in enrichment_results")
            # create color palette if hue values are specified
            if hue_order is not None and hue_values is not None:
                if len(hue_order) != len(hue_values):
                    raise ValueError("Length of hue_order must match length of hue_values")
                palette = dict(zip(hue_order, hue_values))
            elif is_numeric_dtype(top_data[hue_column]):
                nominal_color = len(top_data[hue_column].unique()) < 2

    if top_data.empty:
        return None, None

    # add a newline to the beginning of the size col to provide
    # some separation between the hue and size legends
    nominal_size = True
    if size_column is not None and size_column in top_data.columns:
        new_size_col = "\\n" + size_column
        top_data = top_data.rename(columns={size_column: new_size_col})
        size_column = new_size_col
        if len(top_data[size_column].unique()) > 1:
            nominal_size = False

    # create plots
    fig, ax = plt.subplots()

    # for barplots, use seaborn barplot.
    if plot_type == "barplot":
        # build kwargs to avoid passing color when hue/group is used
        bar_kwargs = dict(data=top_data, y=y_column, x=x_column, orient="h",
                          hue=hue_column, hue_order=hue_order, palette=palette,
                          dodge=False, ax=ax)
        if single_color is not None and hue_column is None:
            bar_kwargs["color"] = single_color

        p = sns.barplot(**bar_kwargs)
        if fig_size is None:
            fig.set_size_inches(5, min(15, top_data.shape[0] / 5))
        ax.set_ylabel("")
        ax.set_title(plot_title)

        try:
            # move legend using the axes object
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        except Exception:
            pass

        return fig, ax

    elif plot_type == "dotplot":
        theme_dict = {**axes_style("whitegrid")}
        theme_dict['axes.prop_cycle'] = cycler(color=["k"])
        p = so.Plot(top_data, x=x_column, y=y_column,
                    color=hue_column, pointsize=size_column)
        p = p.theme(theme_dict)
        p = p.add(so.Dot())

        if nominal_color:
            p = p.scale(color=so.Nominal(values=palette, order=hue_order))
        else:
            p = p.scale(color=so.Continuous(values=palette))

        if nominal_size:
            p = p.scale(pointsize=so.Nominal())
        else:
            p = p.scale(pointsize=so.Continuous(values=dot_sizes))

        p.on(ax).plot()

        ax.grid(visible=False, axis="x")
        ax.grid(visible=True, axis="y")
        ax.set_title(plot_title)
        ax.set_ymargin(1/top_data.shape[0])

        if fig_size is None:
            fig.set_size_inches(5, min(15, top_data.shape[0] / 3))

        if fig.legends:
            fig.legends[0].set_loc("upper left")
            fig.legends[0].set_bbox_to_anchor(
                (1, 1), transform=ax.transAxes)

        return fig, ax

    else:
        raise ValueError(f"Unknown plot type: {plot_type}. Supported plot types are: ['barplot', 'dotplot'].")

# load enrichment results from file
enrichment_results = pd.read_csv(enrichment_file)

if query_name:
    title_prefix = query_name
    prefix = query_name + "_"
else:
    title_prefix = ""
    prefix = ""

if split_columns:
    grouped = enrichment_results.groupby(split_columns, as_index=False)
    for group_vals, group_data in grouped:
        plot_title = [title_prefix] + [": ".join(gv) for gv in zip(split_columns, group_vals)]
        plot_title = "\\n".join(plot_title)
        split_vals = '_'.join([v.strip() for v in group_vals])
        file_prefix = f"{prefix}{split_vals}"
        fig, ax = enrichment_plot(
            group_data, x_column, y_column, plot_type,
            pvalue_column=pvalue_column, pvalue_cutoff=pvalue_cutoff,
            size_column=size_column, dot_sizes=dot_sizes,
            group_column=group_column, group_order=group_order, group_colors=group_colors,
            hue_column=hue_column, hue_order=hue_order, hue_values=hue_values,
            palette=palette, single_color=single_color, sort_ascending=sort_ascending,
            top_n=top_n, plot_title=plot_title, fig_size=fig_size)
        if isinstance(fig, list):
            raise
        if fig is not None:
            fig.savefig(f"{file_prefix}_enrichment_{plot_type}.pdf", bbox_inches="tight")
else:
    group_data = enrichment_results
    if plot_title_column:
        if len(enrichment_results[plot_title_column].unique()) > 1:
            raise ValueError(
            f"Error: plot_title_column '{plot_title_column}' has multiple unique "
            "values but split_columns is not specified. When split_columns is not "
            "specified, title column must have only one unique value to be used for"
            " plot titles.")
        else:
            plot_title = enrichment_results[plot_title_column].iloc[0]
    else:
        plot_title = title_prefix

    file_prefix = prefix
    fig, ax = enrichment_plot(group_data, x_column, y_column, plot_type,
                            pvalue_column=pvalue_column, pvalue_cutoff=pvalue_cutoff,
                            size_column=size_column, dot_sizes=dot_sizes,
                            group_column=group_column, group_order=group_order, group_colors=group_colors,
                            hue_column=hue_column, hue_order=hue_order, hue_values=hue_values,
                            palette=palette, single_color=single_color, sort_ascending=sort_ascending,
                            top_n=top_n, plot_title=plot_title, fig_size=fig_size)
    if fig is not None:
        fig.savefig(f"{file_prefix}enrichment_{plot_type}.pdf", bbox_inches="tight")
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

subprocess.call(["cp", ".command.sh", prefix + "${task.process}.command.sh"])
"""
