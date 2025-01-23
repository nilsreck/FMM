
import sys
import json
import seaborn as sns
from matplotlib import lines
import matplotlib.pyplot as plt
from string import ascii_lowercase
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter


def Calculate_limits(list_mean, list_std):
    lower_error = []
    upper_error = []

    for i in range(len(list_mean)):
        lower = list_mean[i] - list_std[i]
        upper = list_mean[i] + list_std[i]

        if lower < 0:
            lower = list_mean[i]
        else:
            lower = list_std[i]

        if upper > 1:
            upper = 1 - list_mean[i]
        else:
            upper = list_std[i]

        lower_error.append(lower)
        upper_error.append(upper)

    result = [lower_error, upper_error]

    return result


def Plot_Functional_Organization(
    Cancer_type_list,
    Normal_Tissue_list,
    Cell_type_list,
    dimension_list,
    data_path,
    matrix="PPMI",
    Common=False,
    number_similar=500,
    annotation="GO",
    note="Common_Set",
):
    # Paths:

    path_func = data_path + "/"

    if Common == False:
        label = "Individual"
    else:
        label = "Common"

    # Controlers:

    row_control = 0
    label_count = 0

    # Authomatic Labels:

    plot_labels = []

    for i in range(len(Cancer_type_list * 3)):
        plot_labels.append(ascii_lowercase[i])

    # Prepare the grid:

    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({"font.size": 17})
    plt.rc("font", weight="bold")

    fig, axes = plt.subplots(
        len(Cancer_type_list),
        3,
        figsize=(15, 8),
        gridspec_kw={"width_ratios": [1, 1, 0.8]},
    )

    # Start plotting:

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        # Open the file:

        if annotation == "GO":
            Organization = open(
                f"{path_func}Similarity_{tissue}_{note}_{number_similar}.json",
            )
        else:
            Organization = open(
                f"{path_func}Similarity_{tissue}_{note}_{number_similar}_{annotation}.json",
            )

        Organization = json.load(Organization)

        # Line Plot Control Similar:

        axes[row_control, 0].errorbar(
            dimension_list,
            Organization["Sim_Control"],
            yerr=Calculate_limits(
                Organization["Sim_Control"], Organization["Sim_Control_std"]
            ),
            marker="o",
            capsize=5,
            lw=3,
            color="#76de61",
        )

        axes[row_control, 0].errorbar(
            dimension_list,
            Organization["Diss_Control"],
            yerr=Calculate_limits(
                Organization["Diss_Control"], Organization["Diss_Control_std"]
            ),
            marker="p",
            fmt="--",
            capsize=5,
            lw=3,
            color="#76de61",
        )

        axes[row_control, 0].set_xticks(dimension_list)
        axes[row_control, 0].spines["right"].set_visible(False)
        axes[row_control, 0].spines["top"].set_visible(False)
        axes[row_control, 0].spines["left"].set_linewidth(4)
        axes[row_control, 0].spines["bottom"].set_linewidth(4)
        axes[row_control, 0].spines["left"].set_color("grey")
        axes[row_control, 0].spines["bottom"].set_color("grey")
        axes[row_control, 0].spines["bottom"]
        axes[row_control, 0].spines["bottom"]
        axes[row_control, 0].set_ylabel(" ")
        axes[row_control, 0].set_xlabel(" ")
        axes[row_control, 0].set_title(
            #f"{plot_labels[label_count].capitalize()}",
            f"{cancer}",
            fontweight="bold",
            fontsize=17,
            y=0.9,
        )
        label_count = label_count + 1

        if row_control != len(Cancer_type_list) - 1:
            axes[row_control, 0].xaxis.set_major_formatter(NullFormatter())

        # Line plot with dissimilar:

        axes[row_control, 1].errorbar(
            dimension_list,
            Organization["Sim_Cancer"],
            yerr=Calculate_limits(
                Organization["Sim_Cancer"], Organization["Sim_Cancer_std"]
            ),
            marker="o",
            capsize=5,
            lw=3,
            color="#f1948a",
        )
        axes[row_control, 1].errorbar(
            dimension_list,
            Organization["Diss_Cancer"],
            yerr=Calculate_limits(
                Organization["Diss_Cancer"], Organization["Diss_Cancer_std"]
            ),
            marker="p",
            fmt="--",
            capsize=5,
            lw=3,
            color="#f1948a",
        )

        axes[row_control, 1].set_xticks(dimension_list)
        axes[row_control, 1].spines["right"].set_visible(False)
        axes[row_control, 1].spines["top"].set_visible(False)
        axes[row_control, 1].spines["left"].set_linewidth(4)
        axes[row_control, 1].spines["bottom"].set_linewidth(4)
        axes[row_control, 1].spines["left"].set_color("grey")
        axes[row_control, 1].spines["bottom"].set_color("grey")
        axes[row_control, 1].spines["bottom"]
        axes[row_control, 1].spines["bottom"]
        axes[row_control, 1].set_ylabel(" ")
        axes[row_control, 1].set_xlabel(" ")
        axes[row_control, 1].yaxis.set_major_formatter(NullFormatter())
        axes[row_control, 1].set_title(
            #f"{plot_labels[label_count].capitalize()}",
            f"{tissue} control",
            fontweight="bold",
            fontsize=17,
            y=0.9,
        )
        label_count = label_count + 1

        if row_control != len(Cancer_type_list) - 1:
            axes[row_control, 1].xaxis.set_major_formatter(NullFormatter())

        # Bar plot for the tope500 Jaccard:

        ax = sns.barplot(
            x=dimension_list,
            y=Organization["Jaccard_Similar"],
            ax=axes[row_control, 2],
            color="#55a9f3",
        )
        ax.set_yticks([0, 0.95])
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_linewidth(4)
        ax.spines["bottom"].set_linewidth(4)
        ax.spines["left"].set_color("grey")
        ax.spines["bottom"].set_color("grey")
        ax.yaxis.set_major_formatter(NullFormatter())
        if row_control != len(Cancer_type_list) - 1:
            ax.xaxis.set_major_formatter(NullFormatter())
        ax.set_title(
            #f"{plot_labels[label_count].capitalize()}",
            f"{tissue}",
            fontweight="bold",
            fontsize=17,
            y=0.9,
        )

        label_count = label_count + 1

        # Next row:

        row_control = row_control + 1

    fig.text(
        0.55,
        0.05,
        "Dimensions",
        ha="center",
        va="center",
        rotation="horizontal",
        fontsize=20,
        fontweight="bold",
    )

    fig.text(
        0.09,
        0.5,
        "Semantic Similarity",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=20,
        fontweight="bold",
    )

    fig.text(
        0.68,
        0.5,
        "Jaccard Index",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=20,
        fontweight="bold",
    )

    # Generate a manual legend:

    Control_patch = mpatches.Patch(color="#76de61", label="Control")
    Cancer_patch = mpatches.Patch(color="#f1948a", label="Cancer")

    plt.legend(
        handles=[
            Control_patch,
            Cancer_patch,
            lines.Line2D([0], [0], marker="", ls="-", c="black", label="Similar"),
            lines.Line2D([0], [0], marker="", ls="--", c="black", label="Dissimilar"),
        ],
        borderaxespad=0.1,
        bbox_to_anchor=(-0.80, -0.72),
        loc="upper center",
        frameon=True,
        ncol=2,
    )

    if annotation == "GO":
        fig.savefig(
            f"{path_func}Functional_Organization_{note}.svg",
            format="svg",
            dpi=600,
            bbox_inches="tight",
        )
    else:
        fig.savefig(
            f"{path_func}Functional_Organization_{annotation}_{note}.svg",
            format="svg",
            dpi=600,
            bbox_inches="tight",
        )


Plot_Functional_Organization(
    snakemake.params.Cancer_list,
    snakemake.params.Normal_tissue_list,
    snakemake.params.Control_list,
    snakemake.params.dimension_list,
    snakemake.params.data_path,
    matrix="PPMI",
    Common=snakemake.params.use_common,
    number_similar=snakemake.params.number_similar,
    annotation=snakemake.params.annotation,
    note="Common_Set")

