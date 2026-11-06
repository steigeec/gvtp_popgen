#!/usr/bin/env python3
'''
Small-scale PC1 chromplot: no top variance panel, bolder lines, taller
aspect ratio than the stock `winpca chromplot`. Reuses WPCAData/Plot
machinery from the winpca package for data loading, metadata annotation
and color assignment.
'''

import argparse
import sys

WINPCA_DIR = '/private/groups/meyerlab/esteiger/software/winpca'
sys.path.insert(0, WINPCA_DIR)

from modules.data import WPCAData
from modules.plot import Plot
import plotly.graph_objects as go


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('prefix', help='WinPCA run prefix (as used for `winpca pca`)')
    p.add_argument('region', help='chrom:start-end')
    p.add_argument('-p', '--plot_var', default='pc_1', choices=['pc_1', 'pc_2', 'hetp'])
    p.add_argument('-m', '--metadata', required=True)
    p.add_argument('-g', '--groups', required=True)
    p.add_argument('-i', '--interval', type=int, default=None)
    p.add_argument('--numeric', action='store_true')
    p.add_argument('--reverse', action='store_true')
    p.add_argument('-w', '--width', type=int, default=720)
    p.add_argument('-t', '--height', type=int, default=400)
    p.add_argument('--line_width', type=float, default=2.5)
    p.add_argument('--line_opacity', type=float, default=0.55)
    p.add_argument('-o', '--out', required=True, help='output image path (.png/.svg/.pdf)')
    p.add_argument('--legend_out', default=None,
                    help='if set, also write a standalone legend/colorbar image here')
    return p.parse_args()


def main():
    args = parse_args()
    chrom, region = args.region.split(':')
    start, end = (int(x) for x in region.split('-'))

    data = WPCAData(args.prefix)

    plot = Plot(
        args.plot_var,
        prefix=args.prefix,
        data=data,
        chrom=chrom,
        start=start,
        end=end,
        metadata_path=args.metadata,
        color_by=args.groups,
        hex_code_dct=None,
        interval=args.interval,
        numeric=args.numeric,
        reverse=args.reverse,
        plot_fmt_lst=['png'],
    )

    # LOAD & PREPARE DATA (reuse Plot's own logic)
    plot.data_df = getattr(plot.data, plot.plot_var)
    if plot.interval:
        plot.data_df = plot.subset(plot.data_df, plot.interval)
    plot.annotate()
    plot.set_colors()

    # SINGLE-PANEL FIGURE (no top variance-explained strip)
    fig = go.Figure()

    show_legend = not plot.numeric

    for group in reversed(plot.group_lst):
        group_df = plot.data_df[plot.data_df[plot.group_id] == group]

        x_val_lst, y_val_lst = [], []
        for sample in set(group_df['id']):
            sample_df = group_df[group_df['id'] == sample]
            x_val_lst += sample_df['pos'].tolist() + [None]
            y_val_lst += sample_df[plot.plot_var].tolist() + [None]

        fig.add_trace(
            go.Scatter(
                x=x_val_lst,
                y=y_val_lst,
                mode='lines',
                name=str(group),
                legendgroup=str(group),
                line=dict(color=plot.color_dct[group], width=args.line_width),
                opacity=args.line_opacity,
                connectgaps=False,
                showlegend=show_legend,
            )
        )

    # no legend/colorbar in the main plot -- plotted separately (--legend_out)
    fig.update_layout(
        template='simple_white',
        paper_bgcolor='white',
        plot_bgcolor='white',
        font_family='Arial',
        font_color='black',
        width=args.width,
        height=args.height,
        showlegend=False,
        margin=dict(l=80, r=20, t=15, b=40),
    )

    fig.update_xaxes(
        range=[start, end],
        linewidth=1,
        side='bottom', mirror=True,
        ticks='outside', tickfont=dict(size=10), tickformat=',.0f',
        title=None,
    )
    fig.update_yaxes(
        linewidth=1,
        side='left', mirror=True,
        ticks='outside', tickfont=dict(size=16),
        title_font=dict(size=18),
        title=dict(text='<b>' + plot.plot_var_disp, standoff=10),
    )

    fig.write_image(args.out, width=args.width, height=args.height, scale=2)
    print(f'Wrote {args.out}')

    if args.legend_out:
        write_legend(plot, args)


def write_legend(plot, args):
    '''
    Write a standalone legend (categorical) or colorbar (numeric) image,
    decoupled from the main plot so the main plot can use its full width.
    '''

    legend_fig = go.Figure()

    if plot.numeric:
        val_lst = [
            float(x) for x in plot.data_df[plot.color_by]
            if x not in plot.na_lst
        ]
        min_val, max_val = min(val_lst), max(val_lst)
        legend_fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(
                color=[0], coloraxis='coloraxis',
            ),
            showlegend=False,
        ))
        legend_fig.update_layout(
            coloraxis=dict(
                colorscale=plot.color_scale,
                cmin=min_val,
                cmax=max_val,
                colorbar=dict(
                    len=0.9,
                    thickness=20,
                    title=dict(text=plot.color_by, font=dict(size=13), side='right'),
                    tickvals=[min_val, max_val],
                    tickfont=dict(size=12),
                ),
            ),
        )
    else:
        for group in plot.group_lst:
            legend_fig.add_trace(go.Scatter(
                x=[None], y=[None],
                mode='lines',
                name=str(group),
                line=dict(color=plot.color_dct[group], width=args.line_width),
            ))
        legend_fig.update_layout(
            legend=dict(font=dict(size=13)),
        )

    legend_fig.update_layout(
        template='simple_white',
        paper_bgcolor='white',
        plot_bgcolor='white',
        font_family='Arial',
        font_color='black',
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=10, r=10, t=10, b=10),
    )

    legend_fig.write_image(args.legend_out, width=250, height=args.height, scale=2)
    print(f'Wrote {args.legend_out}')


if __name__ == '__main__':
    main()
