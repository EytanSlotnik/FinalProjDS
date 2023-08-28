import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from FinalProjDS.plots23.spliters import apply_split_functions

BP_RANGES = ((0, 49), (50, 59), (60, 64), (65, 69), (70, 74), (75, 79),
             (80, 89), (90, 200))
X_BIN_SIZE = 1
Y_BIN_SIZE = 1
Y_START = {1: 0.5, 2: 0.85, 8: 0.96, 16: 0.985}


# helper function to create plot of nor vs map as heatmap and as scatter plot
# for peaks
def plot_nor_vs_map(df, x_bin_size, y_bin_size, j, rows, y_start):
    map_ = df['cur_bp']
    rate = df['drugrate']
    length = 0.7 / rows
    y_h = y_start - j / (rows - 0.5)
    # Compute the 2D histogram using numpy.histogram2d
    hist, x_edges, y_edges = np.histogram2d(
        map_, rate,
        bins=[
            np.arange(min(map_), 80 + x_bin_size, x_bin_size),
            np.arange(min(rate), max(rate) + y_bin_size, y_bin_size)
        ])

    # Create the go.Histogram2d object
    histogram = go.Histogram2d(
        x=map_,
        y=rate,
        xbins=dict(start=min(map_), end=max(map_), size=x_bin_size),
        ybins=dict(start=min(rate), end=max(rate), size=y_bin_size),
        z=hist,  # Assign the histogram counts to the 'z' property
        colorscale='Viridis',
        colorbar=dict(thickness=15, x=-0.1, y=y_h, len=length))

    # Find and plot max map bin for each nor rate
    max_x_bins = np.argmax(hist, axis=0)
    nor_rate_list = []
    map_list = []
    color_list = []
    sum_bin = hist.sum(axis=0)
    sum_by_rate = np.sum(hist, axis=0)
    max_sum = 0
    # find the maximum x bin for each y bin
    for i, max_x_bin in enumerate(max_x_bins):
        if sum_bin[i] == 0:
            continue
        nor_rate = np.round(y_edges[i], 2)
        map_val = np.round(x_edges[max_x_bin], 2)
        color = hist[max_x_bin][i] * sum_by_rate[i]
        max_sum += color
        nor_rate_list.append(nor_rate)
        map_list.append(map_val)
        color_list.append(color)

    color_list = np.round((color_list / max_sum) * 100, 2)
    hover_text = [f'X: {x_val}<br>Y: {y_val}<br>z: {color_val:.2f}%' for
                  x_val, y_val, color_val in
                  zip(map_list, nor_rate_list, color_list)]

    peak_scatter = go.Scatter(x=map_list,
                              y=nor_rate_list,
                              mode='markers',
                              marker=dict(color=color_list,
                                          colorscale='YlGnBu',
                                          colorbar=dict(title='Percentage',
                                                        thickness=15,
                                                        len=length,
                                                        x=1.1,
                                                        y=y_h
                                                        )),
                              hovertext=hover_text,
                              hoverinfo='text')

    return histogram, peak_scatter


# create plot of nor vs map as heatmap and as scatter plot for peaks for each
# subgroup of the data
def heatmap_and_peak_scatter(df, file_name, split_funcs,
                             title='', x_bin_size=X_BIN_SIZE,
                             y_bin_size=Y_BIN_SIZE):
    dfs, titles = apply_split_functions(df, split_funcs)
    n = len(titles)
    col1_title = 'Heatmap of MAP VS NOR Rate'
    col2_title = 'MAP Peak Per NOR Rate'
    subplot_titles = [col1_title, col2_title] * n

    fig = make_subplots(n, 2, subplot_titles=subplot_titles,
                        row_titles=titles)
    fig.update_layout(height=400 * n, width=1000)

    for i in range(len(titles)):
        heatmap, scatter = plot_nor_vs_map(dfs[i], x_bin_size, y_bin_size, i,
                                           n, Y_START[n])
        fig.add_trace(heatmap, row=i + 1, col=1)
        fig.add_trace(scatter, row=i + 1, col=2)

    fig.update_xaxes(title_text='MAP')
    fig.update_xaxes(tickmode='linear', dtick=1, range=[55, 80], col=2)
    fig.update_yaxes(title_text='NOR RATE', tickmode='linear', dtick=5, col=2)
    fig.update_layout(showlegend=False,
                      title_text=f'Heatmap of MAP VS NOR RATE With Scatter '
                                 f'of Peaks {title}')
    # fig.update_layout(margin=dict(l=100, r=1000, t=100, b=20))

    fig.write_html(f'graphs/heatmaps/{file_name}.html')


# Create bar plot of binned MAP with box plot of drugrate for each bin
def box_plot(data, file_name, split_funcs, title):
    dfs, titles = apply_split_functions(data, split_funcs)
    fig = make_subplots(rows=2, cols=1, vertical_spacing=0.02)

    fig.add_trace(go.Bar(x=titles,
                         y=[len(df) for df in dfs],
                         name='Bin Size'), row=1, col=1)

    # Add box plots for drug rate in each group
    for i, df in enumerate(dfs):
        fig.add_trace(go.Box(y=df['drugrate'],
                             name=f'Drug Rate for MAP {titles[i]}'), row=2,
                      col=1)

    fig.update_layout(xaxis=dict(side='top'),
                      title=title)

    fig.write_html(f'graphs/boxPlots/{file_name}.html')


def corr_plot(data, file_name, split_funcs, title):
    # Create correlation heatmap of MAP and Drugrate for each sub df
    dfs, titles = apply_split_functions(data, split_funcs)
    rows = int(np.ceil(len(titles) / 2))
    fig = make_subplots(rows, 2,
                        subplot_titles=[f'Correlation for {m}' for m in titles]
                        )
    for i, df in enumerate(dfs):
        df = df[['cur_bp', 'drugrate']]

        correlation_matrix = df.corr()
        show = False
        if i == 0:
            show = True
        fig.add_trace(go.Heatmap(z=correlation_matrix.values,
                                 x=correlation_matrix.columns,
                                 y=correlation_matrix.columns,
                                 colorscale='Viridis',
                                 showscale=show),
                      row=i // 2 + 1, col=i % 2 + 1)

    fig.update_layout(title=title)
    fig.write_html(f'graphs/correlation-heatmaps/{file_name}.html')


# Create scatter plot of data with fitted polynomial of selected degrees
def poly_fit_plot(data, file_name, split_funcs, degrees):
    dfs, titles = apply_split_functions(data, split_funcs)
    rows = len(titles)
    mean_fig = make_subplots(rows, 1, subplot_titles=titles)
    var_fig = make_subplots(rows, 1, subplot_titles=titles)
    full_fig = go.Figure()
    mean_count = 0
    var_count = 0
    for df in dfs:
        type_ = df.iloc[0]['type']
        X = df['cur_bp'].values.astype(np.float64)
        y = df['drugrate'].values.astype(np.float64)
        # Create a scatter plot for the original data
        data_trace = go.Scatter(x=X, y=y,
                                mode='markers', name='data')
        if type_ == 'mean':
            data_trace.name = titles[mean_count]
            mean_count += 1
            mean_fig.add_trace(data_trace, row=mean_count, col=1)

        if type_ == 'var':
            data_trace.name = titles[var_count]
            var_count += 1
            var_fig.add_trace(data_trace, row=var_count, col=1)

        full_fig.add_trace(data_trace)

        if df.shape[0] == 0:
            continue
        for degree in degrees:
            print(degree)
            # Fit the polynomial
            coefficients = np.polyfit(X, y, degree)
            coeff_str = ""
            for k, coeff in enumerate(coefficients):
                coeff_str += f"{coeff:.4f}x^{degree - k}"
                if k < degree:
                    coeff_str += " + "

            # Generate fitted data
            fitted_X = np.linspace(min(X), max(X), 100)
            fitted_y = np.polyval(coefficients, fitted_X)
            # Create a line plot for the fitted curve
            fit_trace = go.Scatter(x=fitted_X, y=fitted_y, mode='lines',
                                   name=f'Polynomial Degree '
                                        f'{degree}: {coeff_str}')
            if type_ == 'mean':
                mean_fig.add_trace(fit_trace, row=mean_count, col=1)

            if type_ == 'var':
                var_fig.add_trace(fit_trace, row=var_count, col=1)

            full_fig.add_trace(fit_trace)

    # Update layout
    full_fig.update_yaxes(title_text='Drug Rate')
    mean_fig.update_yaxes(title_text='Drug Rate')
    var_fig.update_yaxes(title_text='Drug Rate')
    full_fig.update_xaxes(title_text='MAP')
    mean_fig.update_xaxes(title_text='MAP')
    var_fig.update_xaxes(title_text='MAP')

    full_fig.update_layout(
        height=300 * rows, width=2000,
        title="Data with Fitted Curve",
        showlegend=True
    )
    mean_fig.update_layout(
        height=300 * rows, width=1500,
        title="Mean of Data with Fitted Curve",
        showlegend=True
    )
    var_fig.update_layout(
        height=3000, width=1500,
        title="Var Of Data with Fitted Curve",
        showlegend=True
    )

    # Show the plot
    full_fig.write_html(f'graphs/fitted-curves/full_{file_name}.html')
    mean_fig.write_html(f'graphs/fitted-curves/mean_{file_name}.html')
    var_fig.write_html(f'graphs/fitted-curves/var_{file_name}.html')
