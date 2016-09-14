#!/usr/bin/env python3.5
# Analyze and produce graphs for TCUP results on mixed cultures
# Fredrik Boulund 2016


from sys import argv, exit
import os
import argparse
import logging

import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt


os.environ["KMP_AFFINITY"] = "disabled"
plt.style.use("ggplot")
logging.basicConfig(level=logging.DEBUG)
pd.options.display.width = None

############################################
## Define the triplicate samples
############################################
# In-solution
triplicates = {"1:1:1:1": ("QE_160819_71", 
                           "QE_160819_74", 
                           "QE_160819_77"),
               "4-SP:2-MC:2-HI:1-SA": ("QE_160819_62", 
                                       "QE_160819_65", 
                                       "QE_160819_68"),
               "4-SA:2-MC:2-HI:1-SP": ("QE_160819_53", 
                                       "QE_160819_56", 
                                       "QE_160819_59")}
# LPI
#triplicates = {"1:1:1:1": ("QE_160819_50", 
#                           "QE_160819_47", 
#                           "QE_160819_44"),
#               "4-SP:2-MC:2-HI:1-SA": ("QE_160819_38", 
#                                       "QE_160819_35", 
#                                       "QE_160819_32"),
#               "4-SA:2-MC:2-HI:1-SP": ("QE_160819_29", 
#                                       "QE_160819_26", 
#                                       "QE_160819_23")}



def parse_commandline(argv):

    parser = argparse.ArgumentParser(description="Analyze and plot TCUP results on mixed cultures")

    parser.add_argument("-s", "--sample-info", dest="SAMPLE_INFO",
            default="sample_info.txt",
            help="Sample information")
    parser.add_argument('-d', '--dpeps-dir', dest="DPEPS_DIR",
            default='samples/',
            help="Folder with discriminative peptides output from TCUP, <sample_name>.discriminative_peptides.txt")
    parser.add_argument('-c', '--correction-factors', dest="CORRECTION_FACTORS",
            default='correction_factors.tab',
            help="Species correction factors")
    parser.add_argument("-m", "--multibars", action="store_true", dest="multibars",
            default=False,
            help="Print multiple bars instead of averages with SEM")
    parser.add_argument("-p", "--show-plots", dest="show_plots", action="store_true",
            default=False,
            help="Show plots using matplotlib.pyplot.show()")

    #if len(argv) < 2:
    #    parser.print_help()
    #    exit(1)
    
    return parser.parse_args()


def parse_sample_info(fn):
    columns = ["Sample", "Treatment", "Species", "Ratio"]
    with open(fn) as f:
        for line in f:
            sample, treatment, mixture, ratios = line.strip().split("\t")
            ratio = ratios.split(":")
            for idx, species in enumerate(mixture.split("-")):
                yield sample, treatment, species, ratio[idx]


def parse_dpeps(fn):
    with open(fn) as f:
        f.readline() # Skip first (sample) line
        f.readline() # Skip second (title) line
        f.readline() # Skip third (header) line
        for line in f:
            cumcount = int(line[:11].strip())
            count = int(line[11:19].strip())
            percentage = float(line[19:25].strip())
            rank = line[25:45].strip()
            description = line[46:].strip()
            if rank == 'species':
                speciesname = " ".join(description.split(maxsplit=2)[:2])
                yield speciesname, cumcount


def parse_dpeps_dir(dpeps_dir, samples):
    for fn in os.scandir(dpeps_dir):
        if fn.is_file() and fn.name.startswith("QE"):
            samplename, content, extension = fn.name.split(".")
            if samplename in samples and content.startswith("taxonomic"):
                fpath = os.path.join(dpeps_dir, fn.name)
                for items in parse_dpeps(fpath):
                    yield (samplename, *items)


def parse_correction_factors(fn):
    correction_factors = {}
    with open(fn) as f:
        for line in f:
            species, factor = line.split("\t")
            correction_factors[species] = float(factor)
    return correction_factors


def analyze_triplicate(mixture_ratio, samples, dpeps, correction_factors, show_plots):
    species = samples.Species.unique()
    sample_colors = samples.Sample.unique()

    # Compute the uncorrected compositions for each replicate
    total_dpeps = dpeps.groupby("Sample").sum()
    compositions = []
    for row in dpeps.itertuples():
        compositions.append(row.Dpeps / total_dpeps.loc[row.Sample])
    dpeps.loc[:, "Proportion"] = compositions

    average_uncorrected_compositions = dpeps.groupby("Species").agg([np.mean, scipy.stats.sem])

    print(mixture_ratio)
    #print(dpeps)
    #print(average_uncorrected_compositions)


    if options.multibars:
        make_factorplot(dpeps,
                "Species",
                "Proportion",
                "Sample",
                "Proportion",
                "Uncorrected proportion",
                "Uncorrected proportions\n"+mixture_ratio,
                mixture_ratio+"_uncorrected")
    else:
        make_barplot(data=average_uncorrected_compositions, 
                column="Proportion",
                title="Uncorrected composition estimation\n"+mixture_ratio,
                filename=mixture_ratio+"_uncorrected",
                )

    
    # Correct the compositions for each species using the correction factors
    new_compositions = []
    for row in dpeps.iterrows():
        new_compositions.append(row[1].Proportion / correction_factors[row[1].Species])
    dpeps.loc[:, "Corrected"] = new_compositions
    #print(dpeps)

    average_corrected_compositions = dpeps.groupby("Species").agg([np.mean, scipy.stats.sem])
    #print(average_corrected_compositions)

    # Normalize within each sample
    sample_normalization = dpeps.groupby("Sample").sum().Corrected
    normalized_proportions = []
    for row, df in dpeps.iterrows():
        normalized_proportions.append(df.Corrected / sample_normalization[df.Sample])
    dpeps.loc[:, "Normalized"] = normalized_proportions
    average_normalized_compositions = dpeps.groupby("Species").agg([np.mean, scipy.stats.sem])
    print(dpeps)
    print(average_normalized_compositions)

    if options.multibars:
        make_factorplot(dpeps, 
                "Species", 
                "Normalized", 
                "Sample", 
                "Proportion", 
                "Normalized proportions", 
                "Normalized estimated composition\n"+mixture_ratio,
                mixture_ratio+"_normalized")
    else:
        make_barplot(data=average_normalized_compositions, 
                column="Normalized",
                title="Normalized composition estimation\n"+mixture_ratio,
                filename=mixture_ratio+"_normalized",
                )
   
    if show_plots:
       plt.show()



def make_barplot(data, column, title, filename):
    fig, ax = plt.subplots(1,1)
    data[column]['mean'].sort_values(ascending=False).plot(kind="bar",
            ax=ax,
            title=title,
            yerr=data[column]['sem'],
            rot=0,
            width=.40,
            color="steelblue",
            error_kw={"elinewidth":1.5,
                      "capsize":4},
            )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_axis_bgcolor("white")
    ax.set_ylim([0,.60])
    #fig.tight_layout()
    fig.savefig(filename+".png")
    fig.savefig(filename+".pdf")
    return fig, ax



def make_factorplot(dpeps, xdata, ydata, huedata, ylabel, xtitle, suptitle, filename):
    import seaborn as sns
    sns.set_style("white")
    sns.set_context("notebook", font_scale=1.0)
    g = sns.factorplot(data=dpeps.sort_values(ydata, ascending=False),
            x=xdata,
            y=ydata,
            hue=huedata,
            kind="bar",
            palette="muted",
            legend=False,
            legend_out=True,
            ci=None, # Disable the bootstrapped CI-based error bar
            #color=(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
            )
    g.despine()
    g.set_ylabels("Proportion")
    g.set_titles("Normalized proportions")
    g.set_xticklabels(rotation=0)
    for axe in g.fig.get_axes():
        new_xticklabels = []
        xticklabels = axe.get_xticklabels()
        for label in xticklabels:
            genus_letter = label.get_text().split()[0][0]+"."
            spname = label.get_text().split()[1]
            new_xticklabels.append(" ".join([genus_letter, spname]))
    g.fig.get_axes()[0].set_xticklabels(new_xticklabels)
    plt.subplots_adjust(top=1.9)
    g.fig.suptitle(suptitle)
    plt.tight_layout()
    g.savefig(filename+".png", 
            dpi=150,
            figsize=[15,10])
    g.savefig(filename+".pdf", 
            dpi=150,
            figsize=[15,10])
    return g




if __name__ == "__main__":
    options = parse_commandline(argv)

    correction_factors = parse_correction_factors(options.CORRECTION_FACTORS)

    sample_columns = ["Sample", "Treatment", "Species", "Ratio"]
    samples = pd.DataFrame(parse_sample_info(options.SAMPLE_INFO), columns=sample_columns)

    dpeps_columns = ["Sample", "Species", "Dpeps"]
    dpeps_generator = parse_dpeps_dir(options.DPEPS_DIR, set(samples.Sample.unique()))
    dpeps = pd.DataFrame(dpeps_generator, columns=dpeps_columns)

    for ratio, triplicate in triplicates.items():
        triplicate_samples = samples[samples.Sample.isin(triplicate)]
        triplicate_dpeps = dpeps[dpeps.Species.isin(triplicate_samples.Species.unique()) & dpeps.Sample.isin(triplicate)].copy()
        analyze_triplicate(ratio, triplicate_samples, triplicate_dpeps, correction_factors, options.show_plots)
