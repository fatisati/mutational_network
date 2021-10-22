import plotly.graph_objects as go


def show_table(df):
    table = go.Table(
        header=dict(values=list(df.index),
                    fill_color='paleturquoise',
                    align='center'),
        cells=dict(values=df.values,
                   fill_color='lavender',
                   align='center'))

    fig = go.Figure(data=[table])
    return fig


def histogram(x, nbins):
    return go.Histogram(x=x, nbinsx=nbins, name=None)


def show_go(fig):
    go.Figure(fig).show()
