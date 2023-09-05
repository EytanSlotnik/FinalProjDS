from FinalProjDS.plots23.spliters import *
import pandas as pd
from FinalProjDS.plots23.creators import *
from sklearn.linear_model import LinearRegression


# Function to remove top 10% of measurements for each patient,
def remove_top_10_percent(group):
    threshold = group['drugrate'].quantile(0.9)
    return group[group['drugrate'] <= threshold]


def filter_drugrate(data, min_, max_):
    return data[(data['drugrate'] < max_) & (data['drugrate'] > min_)]


def filter_map(data, min_, max_):
    return data[(data['cur_bp'] < max_) & (data['cur_bp'] > min_)]


if __name__ == "__main__":
    bp = pd.read_csv('../preprocess/smooth_bp_eicu2.csv')
    max_x = 70
    # Should be part of the general filtering
    bp['drugrate'] = bp['drugrate'].fillna(0)  # fill na with 0
    # remove patients that didn't receive any Nor
    no_nor = (
        bp.groupby('stay_id').apply(lambda x: (x['drugrate']).sum() == 0))
    no_nor = no_nor[no_nor].index.tolist()
    bp = bp[~bp['stay_id'].isin(no_nor)]
    # remove top 10 percent of drugrate per patient
    grouped = bp.groupby('stay_id')
    filtered_bp = grouped.apply(remove_top_10_percent)
    filtered_bp.reset_index(drop=True, inplace=True)
    # filter drugrate
    filtered_drugrate = filter_drugrate(filtered_bp, 1, 50)
    filtered_map = filter_map(filtered_drugrate, 30, max_x)

    # poly_fit_plot(filtered_map, 'binned_poly_fit',
    #               [
    #                   group_by_sections_mean_std],
    #               [2, 3])

    # poly_fit_plot(filtered_map, 'binned_unit_poly_fit',
    #               [
    #                   get_bp_by_unit, group_by_sections_mean_std],
    #               [2, 3])

    poly_fit_plot(filtered_map, f'flipped_unit_change_direction_poly_fit_{max_x}',
                  [
                     get_bp_by_nor_change_with_direction, get_bp_by_unit
                  ],
                  [1, 2, 3], max_x, peak=True)

    # plots:
    # heatmap_and_peak_scatter(filtered_map,
    #                          'MAP_NOR_heatmap_all',
    #                          [default_spliter])
    # heatmap_and_peak_scatter(filtered_map,
    #                          'MAP_NOR_heatmap_unit',
    #                          [get_bp_by_unit],
    #                          title='by Units')
    # heatmap_and_peak_scatter(filtered_map,
    #                          'all_changes_heatmap',
    #                          [get_bp_by_nor_change],
    #                          title='for Change in Rate')
    # heatmap_and_peak_scatter(filtered_map,
    #                          f'heatmap_changes_direction_unit_{max_x}',
    #                          [get_bp_by_nor_change_with_direction,
    #                           get_bp_by_unit],
    #                          max_x,
    #                          title='By Unit For Change In Rate With Direction',
    #                          )
    # box_plot(filtered_drugrate,
    #          'filtered_drug_rate_binned_bar_box',
    #          [get_bp_by_sections],
    #          'Binned MAP With Boxplot for Matching Drugrate (1-50)')
    # box_plot(filtered_bp,
    #          'binned_bar_box',
    #          [get_bp_by_sections],
    #          'Binned MAP With Boxplot for Matching Drugrate')
    # # Correlation map vs nor for all data binned by bp ranges
    # corr_plot(filtered_bp,
    #           'binned_corr_heatmap_1',
    #           [get_bp_by_sections],
    #           'Correlation Heatmap For All data')
    # # # Correlation map vs nor for filtered drugrate 1- 50 binned by bp ranges
    # corr_plot(filtered_drugrate,
    #           'filtered_binned_corr_heatmap_1',
    #           [get_bp_by_sections],
    #           'Correlation Heatmap For Filtered Drugrate (1-50)')
    # # # Correlation Heatmap For All Data Points Where Drugrate Was Changed
    # corr_plot(filtered_bp,
    #           'all_changes_corr_heatmap_with_direction_1',
    #           [get_bp_by_nor_change_with_direction],
    #           'Correlation Heatmap For All Data Points Where Drugrate Was '
    #           'Changed')
    # # # Correlation Heatmap split by units
    # corr_plot(filtered_bp,
    #           'units_corr_heatmap_1',
    #           [get_bp_by_unit],
    #           'Correlation Heatmap For Each Unit')
    # # # Correlation Heatmap for filtered drugrate 1- 50 split by units
    # corr_plot(filtered_bp,
    #           'units_filtered_corr_heatmap_1',
    #           [get_bp_by_unit],
    #           'Correlation Heatmap For Each Unit For Filtered Drugrate (1-50)')
    # Correlation map vs nor for all data binned by bp ranges
    # corr_plot(filtered_bp,
    #           'corr_heatmap_unit_direction',
    #           [get_bp_by_unit, get_bp_by_nor_change_with_direction],
    #           'Correlationn Heatmap By Unit And Change of Nor Rate With Direction')

    # heatmap_and_peak_scatter(filtered_map,
    #                          'MAP_NOR_heatmap_unit_changes_with_direction',
    #                          [get_bp_by_unit,
    #                           get_bp_by_nor_change_with_direction],
    #                          title=
    #                          '<br>by Units For Rows With Changes With Direction')
