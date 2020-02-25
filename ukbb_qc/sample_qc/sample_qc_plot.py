import pandas as pd
from gnomad_hail import *
from gnomad_hail.utils.plotting import *
from statsmodels.robust.scale import mad
from bokeh.layouts import gridplot
from bokeh.plotting import figure
import holoviews as hv

hv.extension("bokeh")

from IPython.display import display_html

relatedness_kin_cutoffs = {"Second-degree": 0.088388, "First-degree": 0.17678}

relatedness_kin_means = {"Second-degree": 0.125, "First-degree": 0.25}


def get_pc_plots(
    pcs_pd,
    pc_name,
    color_col="color",
    colors=None,
    palette=None,
    n_pcs=10,
    plot_height=600,
    plot_width=700,
    alpha=0.7,
) -> Tabs:
    """
    Description
    :param type parameter: Description
    :return: Description of what is returned
    :rtype: Tabs
    """
    plots = []
    factors = pcs_pd[color_col].unique()
    factors.sort()
    if palette is not None:
        palette = palette[len(factors)]
        colors = CategoricalColorMapper(palette=palette, factors=factors)
    for pc in range(1, n_pcs, 2):
        p = figure(title=pc_name)
        legend_it = []
        for factor in factors:
            source = ColumnDataSource(pcs_pd[pcs_pd[color_col] == factor])
            c = p.circle(
                x=f"{pc_name}{pc}",
                y=f"{pc_name}{pc + 1}",
                source=source,
                color={"field": color_col, "transform": colors},
                alpha=alpha,
            )
            legend_it.append((factor, [c]))
        p.xaxis.axis_label = f"PC{pc}"
        p.yaxis.axis_label = f"PC{pc + 1}"
        p.plot_height = plot_height
        p.plot_width = plot_width
        legend = Legend(items=legend_it, location=(10, 10))
        legend.click_policy = "hide"
        p.add_layout(legend, "right")
        plots.append(Panel(child=p, title=f"PC{pc} vs PC{pc + 1}"))

    return Tabs(tabs=plots)


def get_pc_plots_hail(
    pcs_ht,
    pc_name,
    color_col="color",
    colors=None,
    n_pcs=10,
    plot_height=600,
    plot_width=700,
):
    plots = []
    for pc in range(1, n_pcs, 2):
        p = hl.plot.scatter(
            pcs_ht[f"{pc_name}{pc}"],
            pcs_ht[f"{pc_name}{pc + 1}"],
            colors=colors,
            label=pcs_ht[color_col],
            width=plot_width,
            height=plot_height,
            title=pc_name,
            xlabel=f"PC{pc}",
            ylabel=f"PC{pc + 1}",
        )
        plots.append(Panel(child=p, title=f"PC{pc} vs PC{pc + 1}"))
    return Tabs(tabs=plots)


def get_broad_relatedness_plots(
    ibd_pd,
    color_col={"Broad": "Broad_classification"},
    related_factors=["Full-sibling", "Parent-child", "Second-degree", "None"],
    related_colors=["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"],
    ibd_axes=[("0", "1"), ("1", "2"), ("0", "2")],
    data_sources=["Broad"],
):
    colors = CategoricalColorMapper(palette=related_colors, factors=related_factors)

    def plot_ibd(data_source, ibd_1, ibd_2):
        p = figure(title=f"{data_source} IBD{ibd_1} vs IBD{ibd_2}")
        for factor in related_factors:
            source = ColumnDataSource(ibd_pd[ibd_pd[color_col[data_source]] == factor])
            c = p.circle(
                x=f"{data_source}_IBD{ibd_1}",
                y=f"{data_source}_IBD{ibd_2}",
                source=source,
                color={"field": color_col[data_source], "transform": colors},
                alpha=0.6,
                legend=factor,
            )
        p.xaxis.axis_label = f"{data_source} IBD{ibd_1}"
        p.yaxis.axis_label = f"{data_source} IBD{ibd_2}"
        p.legend.click_policy = "hide"
        p.legend.location = "top_right"

        return p

    plots = []
    for ibd_1, ibd_2 in ibd_axes:
        ibd_plots = [
            plot_ibd(data_source, ibd_1, ibd_2) for data_source in data_sources
        ]
        p = gridplot(ibd_plots, ncols=len(ibd_plots), plot_width=500, plot_height=500)
        plots.append(Panel(child=p, title=f"IBD{ibd_1} vs IBD{ibd_2}"))

    return Tabs(tabs=plots)


def get_relatedness_plots(
    ibd_pd,
    color_col={
        "Broad": "Broad_classification",
        "Regeneron": "Regeneron_classification",
    },
    related_factors=["Full-sibling", "Parent-child", "Second-degree", "None"],
    related_colors=["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"],
    ibd_axes=[("0", "1"), ("1", "2"), ("0", "2")],
    data_sources=["Broad", "Regeneron"],
):
    colors = CategoricalColorMapper(palette=related_colors, factors=related_factors)

    def plot_ibd(data_source, ibd_1, ibd_2):
        p = figure(title=f"{data_source} IBD{ibd_1} vs IBD{ibd_2}")
        for factor in related_factors:
            source = ColumnDataSource(ibd_pd[ibd_pd[color_col[data_source]] == factor])
            c = p.circle(
                x=f"{data_source}_IBD{ibd_1}",
                y=f"{data_source}_IBD{ibd_2}",
                source=source,
                color={"field": color_col[data_source], "transform": colors},
                alpha=0.6,
                legend=factor,
            )
        p.xaxis.axis_label = f"{data_source} IBD{ibd_1}"
        p.yaxis.axis_label = f"{data_source} IBD{ibd_2}"
        p.legend.click_policy = "hide"
        p.legend.location = "top_right"

        return p

    plots = []
    for ibd_1, ibd_2 in ibd_axes:
        ibd_plots = [
            plot_ibd(data_source, ibd_1, ibd_2) for data_source in data_sources
        ]
        p = gridplot(ibd_plots, ncols=len(ibd_plots), plot_width=500, plot_height=500)
        plots.append(Panel(child=p, title=f"IBD{ibd_1} vs IBD{ibd_2}"))

    return Tabs(tabs=plots)


def get_outlier_plots(
    outlier_sample_qc,
    meta_ht,
    facet_col="qc_pop",
    qc_metrics=[
        "n_snp",
        "r_ti_tv",
        "r_insertion_deletion",
        "n_insertion",
        "n_deletion",
        "r_het_hom_var",
    ],
):
    # outlier_sample_qc = outlier_sample_qc.annotate(qc_pop=meta_ht[outlier_sample_qc.key].pop.pop,
    #                                               pop_platform_filters=meta_ht[
    #                                                   outlier_sample_qc.key].pop_platform_filters)
    outlier_sample_qc = outlier_sample_qc.annotate(**meta_ht[outlier_sample_qc.key])
    cols = [facet_col]
    key = "s"
    colnames = cols + [key] + [f"sample_qc.{metric}" for metric in qc_metrics]
    new_colnames = cols + [key] + [f"{metric}" for metric in qc_metrics]
    sample_qc_pd = (
        outlier_sample_qc.flatten()
        .select(*colnames)
        .rename(dict(zip(colnames, new_colnames)))
        .key_by(key)
        .to_pandas()
    )
    sample_qc_pd = sample_qc_pd.set_index(key)

    sample_qc_fail_pd = sample_qc_pd.groupby(facet_col).transform(
        lambda x: (x < (x.median() - (4 * mad(x)))) | (x > (x.median() + (4 * mad(x))))
    )
    sample_qc_fail_pd[facet_col] = sample_qc_pd[facet_col]
    plots = None
    tables = []
    factors = sample_qc_pd[facet_col].unique()
    factors.sort()
    for metric in qc_metrics:
        curve_dict = {}
        sample_qc_fail_pd.groupby(facet_col)[metric].value_counts()
        fail_table = (
            sample_qc_fail_pd.groupby(facet_col)[metric]
            .value_counts()
            .unstack()
            .fillna(0)
        )
        fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        tables.append((metric, fail_table))
        for factor in factors:
            source = sample_qc_pd[sample_qc_pd[facet_col] == factor][metric]
            upper = source.median() + (4 * mad(source))
            lower = source.median() - (4 * mad(source))
            p = hv.Distribution(source).opts(
                xrotation=45
            )  # , fill_color={'field': facet_col, 'transform': colors})
            p = p * hv.VLine(lower).opts(line_dash="dashed")
            p = p * hv.VLine(upper).opts(line_dash="dashed")
            curve_dict[factor] = p
        if plots is not None:
            plots = plots + hv.NdLayout(curve_dict, kdims=facet_col, label=metric).opts(
                sizing_mode="scale_both"
            )
        else:
            plots = hv.NdLayout(curve_dict, kdims=facet_col, label=metric).opts(
                sizing_mode="scale_both"
            )

    return tables, plots


def get_multi_outlier_plots(
    outlier_sample_qc,
    meta_ht,
    facet_cols={"pop": "hybrid_pop", "plat": "batch"},
    qc_metrics=[
        "n_snp",
        "r_ti_tv",
        "r_insertion_deletion",
        "n_insertion",
        "n_deletion",
        "r_het_hom_var",
    ],
):

    outlier_sample_qc = outlier_sample_qc.annotate(**meta_ht[outlier_sample_qc.key])
    outlier_sample_qc = outlier_sample_qc.filter(hl.is_defined(outlier_sample_qc.batch))

    cols = list(facet_cols.values())
    key = "s"
    colnames = cols + [key] + [f"sample_qc.{metric}" for metric in qc_metrics]
    new_colnames = cols + [key] + [f"{metric}" for metric in qc_metrics]
    sample_qc_pd = (
        outlier_sample_qc.flatten()
        .select(*colnames)
        .rename(dict(zip(colnames, new_colnames)))
        .key_by(key)
        .to_pandas()
    )
    sample_qc_pd = sample_qc_pd.set_index(key)

    sample_qc_fail_pd = sample_qc_pd.groupby(cols).transform(
        lambda x: (x < (x.median() - (4 * mad(x)))) | (x > (x.median() + (4 * mad(x))))
    )
    sample_qc_fail_pd[cols] = sample_qc_pd[cols]
    plots = None
    tables = []
    pop_factors = sample_qc_pd[facet_cols["pop"]].unique()
    plat_factors = sample_qc_pd[facet_cols["plat"]].unique()
    pop_factors.sort()
    plat_factors.sort()

    for metric in qc_metrics:
        curve_dict = {}
        sample_qc_fail_pd.groupby(["hybrid_pop", "batch"])[
            metric
        ].value_counts().fillna(0)
        fail_table = (
            sample_qc_fail_pd.groupby(cols)[metric].value_counts().unstack().fillna(0)
        )
        # fail_table = fail_table.rename_axis(mapper="None")
        fail_table.columns = ["Pass", "Fail"]
        fail_table["Pct_fail"] = (fail_table["Fail"] / fail_table.sum(axis=1)) * 100
        decimals = pd.Series([0, 0, 2], index=["Pass", "Fail", "Pct_fail"])
        fail_table["Pass"] = pd.to_numeric(fail_table["Pass"], downcast="integer")
        fail_table["Fail"] = pd.to_numeric(fail_table["Fail"], downcast="integer")
        fail_table = fail_table.round(2)  # (decimals)
        fail_table = fail_table.fillna(0)
        tables.append((metric, fail_table))
        for plat in plat_factors:
            for pop in pop_factors:
                factor = f"{plat}_{pop}"
                source = sample_qc_pd[
                    (sample_qc_pd[facet_cols["pop"]] == pop)
                    & (sample_qc_pd[facet_cols["plat"]] == plat)
                ][metric]

                # account for empty groups, otherwise get numpy NaN errors
                if source.empty:
                    continue
                upper = source.median() + (4 * mad(source))
                lower = source.median() - (4 * mad(source))
                p = hv.Distribution(source).opts(xrotation=45)
                p = p * hv.VLine(lower).opts(line_dash="dashed")
                p = p * hv.VLine(upper).opts(line_dash="dashed")
                curve_dict[factor] = p

        if plots is not None:
            plots = plots + hv.NdLayout(
                curve_dict, kdims="batch_population", label=metric
            ).opts(sizing_mode="scale_both")
        else:
            plots = hv.NdLayout(
                curve_dict, kdims="batch_population", label=metric
            ).opts(sizing_mode="scale_both")

    return tables, plots


def kinship_distribution_plot(
    relatedness_ht: hl.Table,
    kin_label: str = "kin",
    degree_cutoffs: Dict[str, float] = relatedness_kin_cutoffs,
    degree_means: Dict[str, float] = relatedness_kin_means,
    plot_height: int = 600,
    plot_width: int = 800,
    label_font_size: int = 16,
    axis_font_size: int = 14,
) -> hv.Histogram:
    """
    Makes a kinship histogram indicating the kinship means and cutoffs for first and
    second degree relatives with vertical lines.

    :param Table relatedness_ht: hail Table containing kinship values between sample pairs
    :param str kin_label: column name containing kinship values
    :param dict of floats degree_cutoffs: Kinship cutoffs for second and first degree relatives
    :param dict of floats degree_means: Kinship means for second and first degree relatives
    :param int plot_height: plot height
    :param int plot_width: plot width
    :param int label_font_size: font size for axis labels and legend text
    :param int axis_font_size: font size for axis tick labels
    :return: a histogram of kinship values
    :rtype: Histogram
    """
    relatedness_pd = relatedness_ht.to_pandas()
    kin = relatedness_pd[kin_label]
    kin = kin[~np.isnan(kin)]
    frequencies, edges = np.histogram(kin, 200)
    p = hv.Histogram((edges, frequencies))
    height = p.data["Frequency"].max() - p.data["Frequency"].min()
    colors = Category10[10]
    for rel, cutoff in degree_cutoffs.items():
        # p = p * hv.VLine(cutoff, label = f'{rel} cutoff').opts(line_dash='dashed')
        p = p * hv.Spikes(
            ([cutoff], [height]), vdims="height", label=f"{rel} cutoff"
        ).opts(line_dash="dashed", color=colors.pop(0), line_width=3)

    for rel, d_mean in degree_means.items():
        # p = p * hv.VLine(d_mean, label = f'{rel} mean')
        p = p * hv.Spikes(
            ([d_mean], [height]), vdims="height", label=f"{rel} mean"
        ).opts(color=colors.pop(0), line_width=3)

    p.opts(
        xlabel="Kinship",
        height=plot_height,
        width=plot_width,
        fontsize={
            "labels": label_font_size,
            "xticks": axis_font_size,
            "yticks": axis_font_size,
        },
    )

    return p


def stacked_bar_AF_proportion(
    ht,
    x_group,
    x_lab,
    y_lab,
    color_group,
    count_group="n",
    count_total="n",
    width=800,
    height=500,
    ylim=(0, 1),
    cmap="Colorblind",
):
    totals_ht = (
        ht.group_by(ht[x_group]).aggregate(bin_total=hl.agg.sum(ht[count_total]))
    ).key_by(x_group)
    frac_ht = ht.group_by(ht[x_group], ht[color_group]).aggregate(
        n_group=hl.agg.sum(ht[count_group])
    )
    frac_ht = frac_ht.annotate(
        fraction=frac_ht.n_group / totals_ht[frac_ht[x_group]].bin_total
    )
    bars = hv.Bars(
        frac_ht.to_pandas(), [(x_group, x_lab), color_group], ("fraction", y_lab)
    )
    bars.opts(
        stacked=True,
        width=width,
        height=height,
        xrotation=45,
        ylim=(0, 1),
        legend_position="right",
        cmap=cmap,
    )

    return bars


def display_tables(table_list):
    head = """<html>
<body>
<table style="width:100%">
<thead align="center">
<tr>
{}
</tr>
</thead>
<tr>
"""
    titles = ""
    row = ""
    for title, serie in table_list:
        s = serie.copy()
        titles += f'<th style="text-align: center;">{title}</th>\n'
        s.name = ""
        row += "<td>{}</td>".format(s.to_html())

    head = head.format(titles)
    head += row
    head += """
</tr>
</table>
</body>
</html>"""
    display_html(head, raw=True)


# NOTE: created for tranche 2, outlier tables were too wide
def display_wide_tables(table_list):
    head = """<html>
<body>
<table style="width:150%">
<thead align="center">
<tr>
{}
</tr>
</thead>
<tr>
"""
    titles = ""
    row = ""
    for title, serie in table_list:
        s = serie.copy()
        titles += f'<th style="text-align: center;">{title}</th>\n'
        s.name = ""
        row += "<td>{}</td>".format(s.to_html())

    head = head.format(titles)
    head += row
    head += """
</tr>
</table>
</body>
</html>"""
    display_html(head, raw=True)
