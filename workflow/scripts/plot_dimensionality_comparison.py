import sys
import json
import matplotlib.pyplot as plt
from string import ascii_lowercase


def Plot_Relative_Error(
    Cancer_type_list,
    Normal_Tissue_list,
    Cell_type_list,
    dimension_list,
    comparison_list,
    data_dir,
    matrix="PPMI",
    annotation="Leaf",
):
    # Path:

    save_cosine = f"{data_dir}/FMM/"
    row_control = 0
    label_count = 0

    # Set the plot characteristics:

    plt.style.use("seaborn-whitegrid")
    plt.rcParams.update({"font.size": 15})
    plt.tight_layout()
    fig, axs = plt.subplots(len(Cancer_type_list), 2, figsize=(15, 10))

    # Reference to the Plots:

    plot_labels = []

    for i in range(len(Cancer_type_list * 2)):
        plot_labels.append(ascii_lowercase[i])

    for cancer, tissue, cell in zip(
        Cancer_type_list, Normal_Tissue_list, Cell_type_list
    ):
        # Load the relatives errors:

        if annotation == "GO":
            with open(
                f"{save_cosine}Relative_Cancer_{cancer}_{matrix}_{annotation}.txt", "r"
            ) as fp:
                Cancer_Error = json.load(fp)
            with open(
                f"{save_cosine}Relative_Control_{tissue}_{cell}_{matrix}_{annotation}.txt",
                "r",
            ) as fp:
                Control_Error = json.load(fp)

        else:
            with open(
                f"{save_cosine}Relative_Cancer_{cancer}_{matrix}_{annotation}.txt", "r"
            ) as fp:
                Cancer_Error = json.load(fp)
            with open(
                f"{save_cosine}Relative_Control_{tissue}_{cell}_{matrix}_{annotation}.txt",
                "r",
            ) as fp:
                Control_Error = json.load(fp)

        # Plot them:

        if len(Cancer_type_list) == 1:
            axs[0].plot(
                comparison_list,
                Cancer_Error,
                marker="o",
                label=f"Cancer {tissue} {matrix}",
                linewidth=5,
                markersize=8,
            )
            axs[0].spines["right"].set_visible(False)
            axs[0].spines["top"].set_visible(False)
            axs[0].spines["left"].set_linewidth(4)
            axs[0].spines["bottom"].set_linewidth(4)
            axs[0].spines["left"].set_color("grey")
            axs[0].spines["bottom"].set_color("grey")
            axs[0].set_ylabel(" ", fontsize=16, fontweight="bold")
            axs[0].set_xlabel(" ", fontsize=16, fontweight="bold")
            axs[0].set_xticks(fontsize=14, rotation=45)
            axs[0].set_yticks(fontsize=14)

            axs[1].plot(
                comparison_list,
                Control_Error,
                marker="o",
                label="Control {Tissue} {matrix}",
                linewidth=5,
                markersize=8,
            )
            axs[1].spines["right"].set_visible(False)
            axs[1].spines["top"].set_visible(False)
            axs[1].spines["left"].set_linewidth(4)
            axs[1].spines["bottom"].set_linewidth(4)
            axs[1].spines["left"].set_color("grey")
            axs[1].spines["bottom"].set_color("grey")
            axs[1].set_ylabel(" ", fontsize=16, fontweight="bold")
            axs[1].set_xlabel(" ", fontsize=16, fontweight="bold")
            axs[1].set_xticks(fontsize=14, rotation=45)
            axs[1].set_yticks(fontsize=14)

        else:
            axs[row_control][0].plot(
                comparison_list,
                Cancer_Error,
                marker="o",
                label=f"Cancer {tissue} {matrix}",
                color="#76de61",
                linewidth=5,
                markersize=8,
            )
            axs[row_control][0].spines["right"].set_visible(False)
            axs[row_control][0].spines["top"].set_visible(False)
            axs[row_control][0].spines["left"].set_linewidth(4)
            axs[row_control][0].spines["bottom"].set_linewidth(4)
            axs[row_control][0].spines["left"].set_color("grey")
            axs[row_control][0].spines["bottom"].set_color("grey")
            axs[row_control][0].set_ylabel(" ", fontsize=16, fontweight="bold")
            axs[row_control][0].set_xlabel(" ", fontsize=16, fontweight="bold")
            axs[row_control][0].set_title(
                plot_labels[label_count].capitalize(),
                fontweight="bold",
                fontsize=17,
                y=1,
            )

            label_count = label_count + 1

            if row_control == 3:
                axs[row_control][0].set_xticklabels(comparison_list, fontsize=14)
                axs[row_control][0].get_xticklabels()[3].set_color("red")
            else:
                axs[row_control][0].set_xticklabels(" ")

            axs[row_control][1].plot(
                comparison_list,
                Control_Error,
                marker="o",
                label="Control {Tissue} {matrix}",
                color="#f1948a",
                linewidth=5,
                markersize=8,
            )
            axs[row_control][1].spines["right"].set_visible(False)
            axs[row_control][1].spines["top"].set_visible(False)
            axs[row_control][1].spines["left"].set_linewidth(4)
            axs[row_control][1].spines["bottom"].set_linewidth(4)
            axs[row_control][1].spines["left"].set_color("grey")
            axs[row_control][1].spines["bottom"].set_color("grey")
            axs[row_control][1].set_ylabel(" ", fontsize=16, fontweight="bold")
            axs[row_control][1].set_xlabel(" ", fontsize=16, fontweight="bold")
            axs[row_control][1].set_title(
                plot_labels[label_count].capitalize(),
                fontweight="bold",
                fontsize=17,
                y=1,
            )

            label_count = label_count + 1

            if row_control == 3:
                axs[row_control][1].set_xticklabels(comparison_list, fontsize=14)
                axs[row_control][1].get_xticklabels()[3].set_color("red")
            else:
                axs[row_control][1].set_xticklabels(" ")

        # Control the rows:

        row_control = row_control + 1

    fig.text(
        0.5,
        0.05,
        "Dimensions",
        ha="center",
        va="center",
        rotation="horizontal",
        fontsize=20,
        fontweight="bold",
    )

    fig.text(
        0.07,
        0.5,
        "RSE",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=20,
        fontweight="bold",
    )

    fig.legend(
        labels=["Control", "Cancer"],
        borderaxespad=0.1,
        bbox_to_anchor=(0.5, 0.01),
        loc="upper center",
        frameon=True,
        ncol=3,
    )

    # Save the plot:

    fig.savefig(f'{data_dir}/Relative_Error_{annotation}.png', format="png", dpi=600,
                 bbox_inches='tight')



# ToDo: Concatenate results from error calculation (optimal dimensionality)

comparison_list = []
dim_list=snakemake.params.dimension_list
for i in range(len(dim_list)-1):
    comparison_list.append(f"{dim_list[i]}-{dim_list[i+1]}")
Plot_Relative_Error(
    snakemake.params.Cancer_list,
    snakemake.params.Normal_tissue_list,
    snakemake.params.Control_list,
    snakemake.params.dimension_list,
    comparison_list,
    snakemake.params.data_path,
    matrix="PPMI",
    annotation=snakemake.params.annotation,
)