import os
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import zscore
from sklearn import decomposition
from sklearn.cluster import KMeans

from pathlib import Path

# Seems to work really well, very neatly divided the GFP+ and GFP- cells into two distinct clusters
# With some expceptions, like 1 cell that is GFP+ but is clustered with GFP- cells
# Oh well, it was never going to be perfect...


def read_data(data_loc, meta_loc, output_loc, plot=True):
    # Some output stuff
    var_output = pd.DataFrame()
    cov_output = pd.DataFrame()

    # Reading data
    passive = pd.read_csv(os.path.join(data_loc, "passive_properties_filtered.csv"))
    spike = pd.read_csv(os.path.join(data_loc, "spike_properties_filtered.csv"))

    # Meaning spike_data
    spike = spike.groupby(["Cell"]).mean(numeric_only=True).reset_index()

    # Merge together spike and passive stuff
    merged_df = pd.merge(passive, spike, on="Cell")

    # Write something to remove excluded cells
    merged_df = merged_df.loc[merged_df.exclude == "N"]

    # Select the columns I want to keep for PCA
    index_cols = ["Cell", "GFP"]
    pca_cols = [
        "Rm",
        "Capacitance",
        "fAHP",
        "Spike_threshold",
        "Spike_amplitude",
        "Half_width",
        "Vm_x",
        "Frequency",
        "Sag ratio",
        "mAHP",
        "Rebound spikes",
        "Firing rate",
    ]
    pca_df = merged_df[pca_cols]
    pca_index = merged_df[
        index_cols
    ].reset_index()  # Will be re-applied after PCA
    # Z-score pca_df
    pca_df = pca_df.apply(zscore)

    # PCA, think this is right, also with indexing and stuff
    pca = decomposition.PCA(n_components=12)
    x = pca.fit_transform(pca_df)
    pca_df = pd.DataFrame(
        {f'{"pca_" + str(i)}': x[:, i] for i in range(pca.n_components_)}
    )
    pca_df["GFP"] = pca_index["GFP"]
    pca_df["Cell"] = pca_index["Cell"]

    # Plot
    if plot:
        fig, ax = plt.subplots()
        for i, tag in enumerate(pca_df["GFP"].unique()):
            d = pca_df.loc[pca_df.GFP == tag]
            if tag == "pos":
                color = "green"
            elif tag == "neg":
                color = "white"
            ax.scatter(
                d["pca_0"],
                d["pca_1"],
                s=20,
                label=tag,
                color=color,
                edgecolors="black",
                linewidths=0.5,
            )
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        fig.tight_layout()
        fig.savefig(os.path.join(output_loc, "PCA_GFP.pdf"))
        plt.show()

    # Explained variance per principal component
    # Explained variance per principal component
    idx = pca.get_feature_names_out()
    var_output.index = idx
    var_output["explained_variance_ratio"] = (
        pca.explained_variance_ratio_ * 100
    )
    fig, ax = plt.subplots(figsize=(2.5, 1.75))
    ax.plot(var_output["explained_variance_ratio"], lw=0.5, color="black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, 60)
    ax.set_yticks([0, 10, 20, 30, 40, 50, 60])
    ax.tick_params(axis="x", labelrotation=0)
    ax.set_xlim(-0.5, 11.5)
    ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    ax.set_xticklabels(
        ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
    )
    ax.set_ylabel("Expl. variance (%)")
    ax.set_xlabel("Principal component")
    ax.grid(lw=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "explvar.pdf"))
    plt.show()

    # Covariance
    cols = pca.feature_names_in_
    cov_output.index = cols
    for i, col in enumerate(cols):
        cov_output[col] = pca.get_covariance()[i]

    # Eigenvalues
    eigenvalues = abs(pca.components_)
    eig_cols = pca.feature_names_in_
    eigen_df = pd.DataFrame(columns=eig_cols)
    for i, v in enumerate(eigenvalues):
        eigen_df.loc[i] = v
    for index, row in eigen_df.iterrows():
        total = sum(row)
        scaled = (row / total) * 100  # Convert to % of total
        # Quick plot of this
        fig, ax = plt.subplots(figsize=(2.5, 1))
        sns.barplot(
            data=scaled.sort_values(),
            color="white",
            edgecolor="black",
            linewidth=0.5,
            width=0.75,
            hatch="////////",
        )

        # Make pretty
        ax.tick_params(axis="x", labelrotation=90)
        ax.xaxis.set_inverted(True)  # or revert sort_values
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylim(0, 15)
        ax.set_ylabel("Variation (%)")
        fig.tight_layout()
        fig.savefig(
            os.path.join(output_loc, "component_var_pca" + str(index) + ".pdf")
        )
        plt.show()

    # Covariance
    cols = pca.feature_names_in_
    cov_output.index = cols
    for i, col in enumerate(cols):
        cov_output[col] = pca.get_covariance()[i]

    # Covariance plots
    if plot:
        fig, ax = plt.subplots()
        corr = cov_output.corr()
        sns.heatmap(corr, ax=ax, vmin=-1, vmax=1, cmap="RdBu")
        fig.get_tight_layout()
        fig.savefig(os.path.join(output_loc, "covariance_plot_gfp.pdf"))
        cov_output.to_excel(
            os.path.join(output_loc, "covariance_matrix_gfp.xlsx")
        )
        var_output.to_excel(
            os.path.join(output_loc, "explained_variance_gfp.xlsx")
        )
    return pca_df


def pca_kmeans(data, output_loc):
    # Outputs
    output_df = pd.DataFrame()
    inertia = []  # Can be used later to find elbow point if u want

    # Some temp indexing, to be re-applied later
    index_cols = ["GFP", "Cell"]
    gfp = data.GFP
    cell = data.Cell
    for c in index_cols:
        data = data.drop(c, axis=1)

    # Make clusters 0-10
    for i in range(10):  # Only 10 clusters, seems reasonable for now
        kmeans = KMeans(n_clusters=i + 1, random_state=0).fit(data)
        for t, label in enumerate(kmeans.labels_):
            output_df.loc[t, "label"] = int(label)
            output_df.loc[t, "pca_1"] = data["pca_0"][t]
            output_df.loc[t, "pca_2"] = data["pca_1"][t]
            output_df.loc[t, "GFP"] = gfp[t]
            output_df.loc[t, "cell"] = cell[t]
        # For now not outputting this to excel

        # Plot for every cluster
        fig, ax = plt.subplots(figsize=(2.5, 2.5))
        # colors = cm.tab20c(np.linspace(0, 1, 2))
        for it, l in enumerate(output_df.label.unique()):
            plot_df = output_df.loc[output_df.label == l]
            plot_df_pos = plot_df.loc[plot_df.GFP == "pos"]
            plot_df_neg = plot_df.loc[plot_df.GFP == "neg"]
            # GFP+ and GFP- cells
            ax.scatter(
                plot_df_pos["pca_1"],
                plot_df_pos["pca_2"],
                label=l,
                lw=0.5,
                color="white",
                hatch="//////////",
                s=15,
            )
            ax.scatter(
                plot_df_neg["pca_1"],
                plot_df_neg["pca_2"],
                lw=0.5,
                color="white",
                s=15,
            )
            # Cluster centers
            centroids = kmeans.cluster_centers_
            for c in range(len(centroids)):
                ax.scatter(
                    centroids[c][0],
                    centroids[c][1],
                    color="red",
                    alpha=0.5,
                    marker="+",
                    s=80,
                )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(-6, 6)
        ax.set_ylim(-4, 4)
        ax.set_xticks([-6, -3, 0, 3, 6])
        ax.set_yticks([-4, -2, 0, 2, 4])
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        # ax.grid(lw=0.25, color="gray")
        fig.tight_layout()
        fig.savefig(os.path.join(output_loc, str(it) + "_cluster_pca.pdf"))
        plt.show()
    return output_df, fig, ax


if __name__ == "__main__":
    # Plotting params
    plt.rcParams["font.sans-serif"] = "Calibri"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.linewidth"] = 0.5
    plt.rcParams["xtick.major.width"] = 0.5
    plt.rcParams["xtick.minor.width"] = 0.5
    plt.rcParams["ytick.major.width"] = 0.5
    plt.rcParams["ytick.minor.width"] = 0.5
    plt.rcParams["hatch.linewidth"] = 0.5
    fonts = {
        "font.size": 8,
        "axes.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.labelsize": 8,
    }
    plt.rcParams.update(fonts)

    # Folders
    # Change data folder to output folder of analysis
    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")

    # Running stuff
    data = read_data(data_loc, meta_loc, output_loc, plot=True)
    pca_kmeans(data, output_loc)
