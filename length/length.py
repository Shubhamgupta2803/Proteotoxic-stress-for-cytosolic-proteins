#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import pingouin as pg

# =========================================================
# COLORS
# =========================================================
BLUE_LIGHT, BLUE_MED, BLUE_DARK = "#84d2f6", "#1f77b4", "#014f86"
ORANGE_LIGHT, ORANGE_MED, ORANGE_DARK = "#ffbf69", "#ff7f0e",  "#cc5803"

def clean_df(df):
    df = df.loc[:, ~df.columns.str.contains("^Unnamed")]
    return df.reset_index(drop=True)

def group_to_label(g):
    thr = int(g["thr"]) if g["thr"].is_integer() else g["thr"]
    return f"{g['kind']}_{g['mito']}_{g['op']}{thr}"

def resolve_file(group, folder, kind):
    thr = int(group["thr"]) if group["thr"].is_integer() else group["thr"]
    name = f"{group['kind']}_{group['mito']}_{group['op']}{thr}_{kind}.csv"
    path = os.path.join(folder, name)
    return path if os.path.exists(path) else None

def extract_op_thr(label):
    label = label.strip()
    parts = label.split("_")
    if len(parts) < 3:
        raise ValueError(f"Invalid label format: {label}")
    opthr = parts[2]
    if opthr.startswith(">="):
        return ">=", float(opthr[2:])
    elif opthr.startswith("<"):
        return "<", float(opthr[1:])
    else:
        raise ValueError(f"Invalid operator in label: {label}")

def build_threshold_map(labels):
    threshold_map = {}
    for label in labels:
        label = label.strip()
        prefix = label[0]
        op, val = extract_op_thr(label)
        if op == ">=":
            threshold_map.setdefault(prefix, []).append(val)
    for prefix in threshold_map:
        threshold_map[prefix] = sorted(set(threshold_map[prefix]))
    return threshold_map

def get_color(label, threshold_map):
    label = label.strip()
    op, val = extract_op_thr(label)
    prefix = label[0]
    if prefix == "I":
        light, med, dark = BLUE_LIGHT, BLUE_MED, BLUE_DARK
    elif prefix == "T":
        light, med, dark = ORANGE_LIGHT, ORANGE_MED, ORANGE_DARK
    else:
        raise ValueError(f"Unknown prefix: {prefix}")
    if op == "<":
        return light
    thresholds = threshold_map.get(prefix, [])
    if not thresholds:
        return med
    if len(thresholds) == 1:
        return med
    min_thr = min(thresholds)
    max_thr = max(thresholds)
    if val == min_thr:
        return med
    elif val == max_thr:
        return dark
    else:
        return med

def plot_ecdf(dfA, dfB, labelA, labelB, metric, out_png, pdf_obj=None, k_val=None, threshold_map=None):
    a_vals = pd.to_numeric(dfA[metric], errors="coerce").dropna()
    b_vals = pd.to_numeric(dfB[metric], errors="coerce").dropna()
    if a_vals.empty:
        a_vals = pd.Series([0])
    if b_vals.empty:
        b_vals = pd.Series([0])

    def pretty_label(raw):
        parts = raw.split("_")
        kind = "Interactome" if parts[0] == "I" else "Translatome"
        return f"{kind} {parts[2]}%"

    groupA = pretty_label(labelA)
    groupB = pretty_label(labelB)

    ecdf_df = pd.concat([
        pd.DataFrame({"Group": groupA, metric: a_vals}),
        pd.DataFrame({"Group": groupB, metric: b_vals})
    ], ignore_index=True)

    colA = get_color(labelA, threshold_map)
    colB = get_color(labelB, threshold_map)
    palette = {groupA: colA, groupB: colB}

    fig, ax = plt.subplots(figsize=(5, 4.2))
    sns.ecdfplot(data=ecdf_df, x=metric, hue="Group", palette=palette, linewidth=1.5, ax=ax)

    if metric == "Length":
        display_metric = "Sequence Length"
    else:
        display_metric = metric.replace("_", " ").title()

    type_str = "Mito" if "_mito_" in labelA.lower() else "Non-Mito"
    try:
        power_val = int(float(k_val) * 100)
    except:
        power_val = k_val

    title = f"{display_metric} ECDF for {type_str} Proteins in $L_{{{power_val}}}$"
    ax.set_title(title, fontsize=10, pad=12)
    ax.set_xlabel(display_metric, fontsize=12)
    ax.set_ylabel("Cumulative Probability", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    sns.move_legend(ax, "center right", bbox_to_anchor=(1.05, 0.5), frameon=False, title=None, fontsize=7)

    plt.tight_layout()
    base = os.path.splitext(out_png)[0]
    ecdf_out = f"{base}_ecdf.svg"
    fig.savefig(ecdf_out, format="svg", bbox_inches="tight")
    if pdf_obj:
        pdf_obj.savefig(fig)
    plt.close(fig)

def plot_metric(dfA, dfB, labelA, labelB, metric, out_png, pdf_obj, k_val, threshold_map):
    if metric not in dfA.columns or metric not in dfB.columns: return
    dfA, dfB = clean_df(dfA[dfA[metric].notna()]), clean_df(dfB[dfB[metric].notna()])
    if dfA.empty or dfB.empty: return

    actual_max = max(dfA[metric].max(), dfB[metric].max())
    use_log = (metric == "Length" or actual_max > 1500)

    try:
        mwu = pg.mwu(dfA[metric], dfB[metric])
        p, r = mwu["p-val"].values[0], mwu["RBC"].values[0]
    except: p, r = 1.0, 0.0

    colA = get_color(labelA, threshold_map)
    colB = get_color(labelB, threshold_map)

    def pretty_label(raw):
        p = raw.split("_")
        kind = 'Interactome' if p[0]=='I' else 'Translatome'
        return f"{kind} {p[2]}%"

    df_plot = pd.DataFrame({
        metric: pd.concat([dfA[metric], dfB[metric]], ignore_index=True),
        "Group": [pretty_label(labelA)]*len(dfA) + [pretty_label(labelB)]*len(dfB)
    })

    fig, ax = plt.subplots(figsize=(4, 4.4))
    if use_log:
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
        ax.yaxis.set_minor_locator(ticker.NullLocator())

    sns.boxplot(data=df_plot, x="Group", y=metric, palette=[colA, colB],
                width=0.7, linewidth=1, ax=ax, showfliers=True,
                flierprops=dict(marker='o', markerfacecolor='none', markeredgecolor='black',
                                markeredgewidth=0.5, alpha=0.8, markersize=4))

    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=10)

    y_min, y_max = ax.get_ylim()
    if use_log:
        ax.set_ylim(bottom=y_min*0.125, top=y_max*6)
    else:
        ax.set_ylim(bottom=y_min - (y_max-y_min)*0.03, top=y_max + (y_max-y_min)*0.125)

    bracket_y = 0.93
    ax.plot([0, 0, 1, 1], [bracket_y-0.02, bracket_y, bracket_y, bracket_y-0.02],
            transform=ax.get_xaxis_transform(), color="black", linewidth=1, clip_on=False)
    if p < 1e-3:
        p_text = f"{p:.3e}"
    else:
        p_text = f"{p:.3f}"

    ax.text(0.5, bracket_y + 0.01, f"p = {p_text}", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=8)

    means = df_plot.groupby("Group", sort=False)[metric].mean()
    for i, m in enumerate(means):
        ax.hlines(m, i-0.2, i+0.2, linestyles="--", colors="black", linewidth=1)

    ax.text(0.02, 0.98, f"-- Mean\nr ={r:.3f}", transform=ax.transAxes, ha="left", va="top", fontsize=9)
    for i, n in enumerate(df_plot.groupby("Group", sort=False).size()):
        ax.text(i, 0.04, f"n = {n}", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8)

    if metric == "Length":
        display_metric = "Sequence Length"
    else:
        display_metric = metric.replace("_", " ").title()

    type_str = "Mito" if "_mito_" in labelA.lower() else "Non-Mito"
    try:
        power_val = int(float(k_val) * 100)
    except:
        power_val = k_val

    full_title = f"{display_metric} for {type_str} Proteins in $L_{{{power_val}}}$"
    ax.set_title(full_title, fontsize=10, pad=12)
    ax.set_ylabel(display_metric, fontsize=12)
    ax.set_xlabel("")

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    fig.savefig(out_png.replace(".png", ".svg"))
    if pdf_obj: pdf_obj.savefig(fig)
    plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotated", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--k_fraction", type=str, default="0.9")
    parser.add_argument("--compare")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    pdf = PdfPages(os.path.join(args.output, "ALL_PLOTS.pdf"))

    tasks = [
        (args.annotated, "annotated", ["Length"]),
    ]

    for folder, kind, metrics in tasks:
        if not os.path.exists(folder): continue

        files = sorted([f for f in os.listdir(folder) if f.startswith("I_") and f.endswith(f"_{kind}.csv")])
        all_labels = []

        for f in files:
            p = f.replace(f"_{kind}.csv", "").split("_")
            gI = {"kind":"I", "mito":p[1], "op":(">=" if ">=" in p[2] else "<"),
                  "thr":float(p[2].replace(">=","").replace("<",""))}
            gT = gI.copy(); gT["kind"] = "T"
            all_labels.append(group_to_label(gI))
            all_labels.append(group_to_label(gT))

        threshold_map = build_threshold_map(all_labels)

        for f in files:
            p = f.replace(f"_{kind}.csv", "").split("_")
            gI = {"kind":"I", "mito":p[1], "op":(">=" if ">=" in p[2] else "<"), "thr":float(p[2].replace(">=","").replace("<",""))}
            gT = gI.copy(); gT["kind"] = "T"
            ip, tp = resolve_file(gI, folder, kind), resolve_file(gT, folder, kind)
            if ip and tp:
                dfI, dfT = pd.read_csv(ip), pd.read_csv(tp)
                for metric in metrics:
                    base_name = f"{metric}_{group_to_label(gI)}_VS_{group_to_label(gT)}"
                    box_out = os.path.join(args.output, f"{base_name}.png")
                    ecdf_out = os.path.join(args.output, f"{base_name}.png")
                    plot_metric(dfI, dfT, group_to_label(gI), group_to_label(gT), metric, box_out, pdf, args.k_fraction, threshold_map)
                    plot_ecdf(dfI, dfT, group_to_label(gI), group_to_label(gT), metric, ecdf_out, pdf, args.k_fraction, threshold_map)

        if args.compare:
            try:
                comparisons = args.compare.split(";")
                for comp in comparisons:
                    lA, lB = comp.split(",")
                    fA = os.path.join(folder, f"{lA}_{kind}.csv")
                    fB = os.path.join(folder, f"{lB}_{kind}.csv")
                    if os.path.exists(fA) and os.path.exists(fB):
                        dfA = pd.read_csv(fA)
                        dfB = pd.read_csv(fB)
                        for metric in metrics:
                            base_name = f"COMPARE_{metric}_{lA}_VS_{lB}"
                            box_out = os.path.join(args.output, f"{base_name}.png")
                            ecdf_out = os.path.join(args.output, f"{base_name}.png")
                            plot_metric(dfA, dfB, lA, lB, metric, box_out, pdf, args.k_fraction, threshold_map)
                            plot_ecdf(dfA, dfB, lA, lB, metric, ecdf_out, pdf, args.k_fraction, threshold_map)
            except Exception as e:
                print(f"Compare Error: {e}")

    pdf.close()