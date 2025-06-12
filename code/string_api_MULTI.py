#!/usr/bin/env python3

"""This script analyses a set of gene/proteins with STRING DB and
gets back a set of plots and tables useful for posterior analyses."""

# Author: Daniel Martinez-Martinez
# Year: 2021 (Pandemic times)

from time import sleep
from collections import Counter

import os
import json
import argparse
import pandas as pd
import requests
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


# define functions

def gene_list(file_list):
    """
    reads a gene list from a txt files
    """
    data = pd.read_csv(file_list, sep=' ', header=None)
    genes = data[0].values.tolist()
    return genes


def get_samples(ex_file):
    """
    For multi-sample files, gets
    """
    sheet_list = list(pd.read_excel(ex_file, None).keys())
    list_of_samples = []
    for name in sheet_list:
        temp = name.split('_')[:-1]
        temp = '_'.join(temp)
        list_of_samples.append(temp)
        list_of_samples = list(set(list_of_samples))
    return list_of_samples


def get_net_image(genes, species=511145, out_net='full_network.svg'):
    """
    This function gets a gene list as an input and
    outputs a svg image of the network from those genes
    from string- db
    A different species and output name can be chosen
    """
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "svg"
    method = "network"

    # Construct URL
    request_url = "/".join([string_api_url, output_format, method])

    # Parameters
    params = {
        "identifiers": "\r".join(genes),  # your protein
        "species": species,  # species NCBI identifier
        "network_flavor": "confidence",  # show confidence links
    }

    response = requests.post(request_url, data=params)

    print(f"Saving interaction network to {out_net}.svg file")

    with open(f'./{sub_folder}/{out_net}.svg', 'wb') as fh_net:
        fh_net.write(response.content)

    sleep(1)


def get_enrichment_data(genes, species=511145):
    """
    Function gets gene list and extracts functional enrichment (if any)
    """
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "json"
    method = "enrichment"

    # Construct the request
    request_url = "/".join([string_api_url, output_format, method])

    # Set parameters
    params = {

        "identifiers": "%0d".join(genes),  # your protein
        "species": species,  # species NCBI identifier
        "caller_identity": "www.awesome_app.org"  # your app name

    }

    # Call STRING
    response = requests.post(request_url, data=params)
    # Read the data
    data = json.loads(response.text)
    # transform data to a dataframe
    data_long = pd.DataFrame(data)
    return data_long


def count_words(df, category='Process', nwords=10):
    """
    This function inputs the enrichment dataframe from get_enrichment_data
    and outputs a list of names with their relative frequencies
    It can be parsed per category type
    The function controls for useless words (e.g. to, and, of...) and for
    other punctuation symbols
    """
    proc_names = df[df['category'] == category]['description'].tolist()
    words = []
    for elm in proc_names:
        for word in elm.split():
            word = word.replace(',', '')
            words.append(word.lower())
    # create a data frame
    word_df = pd.DataFrame.from_dict(Counter(words), orient='index', columns=['count'])
    # filter words and keep the top ones
    filter_list = ['process', 'substance', 'to', 'a', 'metabolic', 'via',
                   'and', 'of', 'incl.', 'by', 'in', 'with', 'the', 'from', 'i',
                   'on', 'for']
    word_df = word_df[~word_df.index.isin(filter_list)]
    word_df = word_df.sort_values('count', ascending=False)
    # calculate the relative word use in this sublist
    word_df = word_df.head(nwords)
    word_df = word_df.assign(relative=word_df['count'] / word_df['count'].sum())
    return word_df


def radar_chart_single(dft, category):
    """
    This function inputs the word count resulting from count_words
    and plots the radar chart from that list
    """
    # transpose the dataframe
    dft = dft.transpose()

    # Each attribute we'll plot in the radar chart.
    labels = list(dft)[0:]
    # Let's take the values
    values = dft.loc['relative'].tolist()

    # Number of variables we're plotting.
    num_vars = len(labels)

    # Split the circle into even parts and save the angles,
    # so we know where to put each axis.
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

    # The plot is a circle, so we need to "complete the loop"
    # and append the start value to the end.
    values += values[:1]
    angles += angles[:1]

    # ax = plt.subplot(polar=True)
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

    # Draw the outline of our data.
    ax.plot(angles, values, color='red', linewidth=1)
    # Fill it in.
    ax.fill(angles, values, color='red', alpha=0.25)

    # Fix axis to go in the right order and start at 12 o'clock.
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    # Draw axis lines for each angle and label.
    ax.set_thetagrids(np.degrees(angles[:-1]), labels)

    # Go through labels and adjust alignment based on where
    # it is in the circle.
    for label, angle in zip(ax.get_xticklabels(), angles):
        if angle in (0, np.pi):
            label.set_horizontalalignment('center')
        elif 0 < angle < np.pi:
            label.set_horizontalalignment('left')
        else:
            label.set_horizontalalignment('right')

    ax.set_rlabel_position(180 / num_vars)  # sets y-label to middle

    # Add some custom styling.
    ax.tick_params(colors='#222222')  # Change the color of the tick labels.
    ax.tick_params(axis='y', labelsize=8)  # Make the y-axis (0-100) labels smaller.
    ax.grid(color='#AAAAAA')  # Change the color of the circular gridlines.
    ax.spines['polar'].set_color('#222222')  # Change color of the outermost gridline
    ax.set_facecolor('#FAFAFA')  # circle background color
    fig.tight_layout()
    # Lastly, give the chart a title and give it some
    # padding above the "Acceleration" label.
    ax.set_title(category, y=1.08)


def get_multi_table(up, down, cat, nwords=10):
    """
    Inputs two enrichment tables from STRING and outputs a
    table with the relative count of the most common words 
    within a specified category.
    """
    # get up and down words
    up_words = count_words(up, category=cat, nwords=nwords)
    down_words = count_words(down, category=cat, nwords=nwords)

    up_words = up_words.drop('count', axis=1)
    up_words.columns = ['up']
    down_words = down_words.drop('count', axis=1)
    down_words.columns = ['down']

    # join both dataframes, fill NA with 0s
    total_words = up_words.join(down_words, how='outer', sort=True).fillna(0)

    # sorting values will make the plot look better
    total_words = total_words.sort_values(['up', 'down'])
    return total_words


def add_to_radar(dft, direction, color, ax, angles):
    """
    Inner function used within radar_chart_multi
    Gets the data in shape to include it in the radar plot
    """
    values = dft.loc[direction].tolist()
    values += values[:1]
    ax.plot(angles, values, color=color, linewidth=1, label=direction)
    ax.fill(angles, values, color=color, alpha=0.25)


def radar_chart_multi(gene_table, category):
    """
    Takes a table with two columns (named up and down) from the function get_multi_table,
    and returns a spider plot from it
    """
    # transpose the dataframe
    dft = gene_table.transpose()

    # Each attribute we'll plot in the radar chart.
    labels = list(dft)[0:]

    # Number of variables we're plotting.
    num_vars = len(labels)

    # Split the circle into even parts and save the angles,
    # so we know where to put each axis.
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

    # The plot is a circle, so we need to "complete the loop"
    # and append the start value to the end.
    angles += angles[:1]

    # ax = plt.subplot(polar=True)
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

    # add samples to the radar plot
    add_to_radar(dft, 'up', '#1aaf6c', ax, angles)
    add_to_radar(dft, 'down', '#429bf4', ax, angles)

    # Fix axis to go in the right order and start at 12 o'clock.
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    # Draw axis lines for each angle and label.
    ax.set_thetagrids(np.degrees(angles[:-1]), labels)

    for label, angle in zip(ax.get_xticklabels(), angles):
        if angle in (0, np.pi):
            label.set_horizontalalignment('center')
        elif 0 < angle < np.pi:
            label.set_horizontalalignment('left')
        else:
            label.set_horizontalalignment('right')

    # set limits
    y_lim = gene_table.max().max() * 1.05
    ax.set_ylim(0, y_lim)
    ax.set_rlabel_position(180 / num_vars)  # set position of y labels

    # Add some custom styling.
    # Change the color of the tick labels.
    ax.tick_params(colors='#222222')
    # Make the y-axis (0-100) labels smaller.
    ax.tick_params(axis='y', labelsize=8)
    # Change the color of the circular gridlines.
    ax.grid(color='#AAAAAA')
    ax.spines['polar'].set_color('#222222')  # Change color of outermost gridline
    ax.set_facecolor('#FAFAFA')
    ax.set_title(category, y=1.08)  # Add title.
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))  # legend
    fig.tight_layout()


# define options to parse
# Create the parser
my_parser = argparse.ArgumentParser(
    prog='STRING API enrich',
    description='List the content of a folder')

# Add the arguments
my_parser.add_argument('Input',
                       metavar='-i',
                       type=str,
                       help='input file')
my_parser.add_argument('Output',
                       metavar='-o',
                       type=str,
                       help='output name for output files')
my_parser.add_argument('Species',
                       metavar='-s',
                       type=str,
                       help='select between ecoli or human')

# Execute the parse_args() method
args = my_parser.parse_args()

# redefine args
filename = args.Input
output = args.Output

# try creation of folder
try:
    os.mkdir(output)
except OSError:
    print(f"Creation of the directory {output} failed")
else:
    print(f"Successfully created the directory {output}")

# program version
_VERSION_ = 0.3

# define the list of species included in the script
species_list = {
    'ecoli': 511145,
    'human': 9606,
    'celegans': 6239,
    'mouse': 10090,
    'fly': 7227,
    'zebrafish': 7955,
    'yeast': 4932,
    'PA': 208964
}

spc = species_list[args.Species]


# MAIN PART
def main():
    """
    Main function of the script
    """
    print(f'\nAnalysing the file {filename}\n')
    # read tables from excel, see how many samples are
    samples = get_samples(filename)

    print(f'The file has these samples{samples}\n')

    # for each sample and for each direction, get the network, enrichment file and enrichment summary plot
    for sample in samples:

        global sub_folder
        sub_folder = output + '/' + sample
        # try creation of folder
        try:
            os.mkdir(sub_folder)
        except OSError:
            print(f"Creation of the directory {sample} failed, probably it already exists?")
        else:
            print(f"Successfully created the directory {sample}")

        up = pd.read_excel(filename, sample + '_UP')
        down = pd.read_excel(filename, sample + '_DOWN')

        # create lists of genes
        up_genes = up.iloc[:, 0].tolist()
        down_genes = down.iloc[:, 0].tolist()

        # get networks (it works fine)
        print('Getting the network for the UPregulated elements\n')
        get_net_image(up_genes, species=spc, out_net=f'{sample}_up_network.svg')
        print('Getting the network for the DOWNregulated elements\n')
        get_net_image(down_genes, species=spc, out_net=f'{sample}_down_network.svg')

        # see how many categories are shared between up and down, plot the radar chart
        # get enrichment
        print(f'Getting enrichment for sample {sample}\n')
        global up_enrich, down_enrich  # define global variables within function
        up_enrich = get_enrichment_data(up_genes, species=spc)
        down_enrich = get_enrichment_data(down_genes, species=spc)


        # test that we have enrichment data, if not, pass
        if up_enrich.shape[0] > 0 or down_enrich.shape[0] > 0:

            # test whether the categories are shared or not
            if up_enrich.shape[0]:
                cat_up = set(up_enrich['category'].unique().tolist())
            else:
                cat_up = set()

            if down_enrich.shape[0]:
                cat_dw = set(down_enrich['category'].unique().tolist())
            else:
                cat_dw = set()

            # all the shared categories
            shared_cats = cat_up.intersection(cat_dw)

            if len(shared_cats) > 0:
                print('\nPlotting shared categories as radar plots!\n')
                print(shared_cats)
                print('\n')

                for category in list(shared_cats):
                    word_table = get_multi_table(up_enrich, down_enrich, cat=category, nwords=10)
                    radar_chart_multi(word_table, category)
                    plt.savefig(f'./{sub_folder}/{sample}_{category}_radar_chart.pdf')

            # if there are categories not present in both, plot separate plots for each of them 
            if cat_up.difference(cat_dw) != set():
                print(f'Single categories {cat_up.difference(cat_dw)} were found for the UP case!\n')
                print('Plotting them!\n')
                for cat in cat_up.difference(cat_dw):
                    word_df = count_words(up_enrich, category=cat, nwords=10)
                    radar_chart_single(word_df, category=cat)
                    plt.savefig(f'./{sub_folder}/{sample}_{cat}_UP_radar_chart.pdf')
            else:
                pass

            if cat_dw.difference(cat_up) != set():
                print(f'Single categories {cat_dw.difference(cat_up)} were found for the DOWN case!\n')
                print('Plotting them!\n')
                for cat in cat_dw.difference(cat_up):
                    word_df = count_words(down_enrich, category=cat, nwords=10)
                    radar_chart_single(word_df, category=cat)
                    plt.savefig(f'./{sub_folder}/{sample}_{cat}_DOWN_radar_chart.pdf')

            else:
                pass
                # print(f'No single category was found for sample {sample}!\n')

            # join both datasets into one and save it

            up_enrich['direction'] = 'UP'
            down_enrich['direction'] = 'DOWN'
            enrich = pd.concat([up_enrich, down_enrich], axis = 0)

            print(f'\nSaving enrichment in file {sample}_output.xlsx\n')
            with pd.ExcelWriter(f'./{sub_folder}/{sample}_output.xlsx') as writer:
                for element in enrich.category.unique():
                    enrich_df = enrich[enrich['category'] == element]
                    enrich_df.to_excel(writer, sheet_name=element)
        elif up_enrich.shape[0] == 0:
            print(f'There was not enrichment for Sample {sample}!!')
            pass

    print('\nAll analyses have finished!\n')


if __name__ == '__main__':
    main()