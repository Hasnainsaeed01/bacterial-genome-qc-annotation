#!/usr/bin/env python3
"""
Genome QC + annotation summary for 10 Salmonella enterica assemblies.

Reads:
  results/quality_report.tsv   CheckM2 predict output (completeness, contamination, ...)
  prokka_out/<genome>/<genome>.txt   Prokka per-genome stats (CDS, tRNA, rRNA, ...)
  data/genome_metadata.csv     accession / serovar / strain labels

Writes:
  results/genome_qc_summary.csv   one row per genome with QC verdict
  figures/01_completeness_contamination.png
  figures/02_features_per_genome.png
  figures/03_size_gc_per_genome.png

QC rule (stated in README.md): a genome is "high quality" when
CheckM2 completeness >= 90% and contamination <= 5%.
"""
import os
import re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(REPO, "results")
FIGURES = os.path.join(REPO, "figures")
os.makedirs(FIGURES, exist_ok=True)

# 1. CheckM2 report. The user ran CheckM2 predict on his own machine
#    (v1.0.1, 8 threads) on the same 10 assemblies; the resulting
#    checkm2_quality_report.tsv is bundled in results/.
qr_path = os.path.join(RESULTS, "checkm2_quality_report.tsv")
if not os.path.exists(qr_path):
    raise SystemExit("results/checkm2_quality_report.tsv not found. "
                     "See README.md for how it was produced.")
checkm = pd.read_csv(qr_path, sep="\t")
checkm = checkm.rename(columns={"Name": "genome"})
# CheckM2 names rows after the input file basename; strip extension just in case.
checkm["genome"] = checkm["genome"].astype(str).str.replace(r"\.fna$", "", regex=True)

# 2. Prokka stats: parse "<key>: <value>" lines from each <genome>.txt.
#    Only the small .txt stats files are bundled (results/prokka_stats/);
#    full Prokka output dirs stay in the workdir (too heavy for GitHub)
#    and are reproduced by run_qc_annotation.sh.
prokka_rows = []
prokka_dir = os.path.join(RESULTS, "prokka_stats")
for fname in sorted(os.listdir(prokka_dir)):
    if not fname.endswith(".txt"):
        continue
    genome = fname[:-4]
    txt = os.path.join(prokka_dir, fname)
    if not os.path.isfile(txt):
        continue
    stats = {"genome": genome}
    with open(txt) as fh:
        for line in fh:
            m = re.match(r"\s*([^:]+):\s*([\d.]+)", line)
            if m:
                key = m.group(1).strip().lower().replace(" ", "_")
                try:
                    stats[key] = int(float(m.group(2)))
                except ValueError:
                    pass
    prokka_rows.append(stats)
prokka = pd.DataFrame(prokka_rows)

# 3. Metadata labels (accession, serovar).
meta = pd.read_csv(os.path.join(REPO, "data", "genome_metadata.csv"))
meta = meta.rename(columns={"label": "genome"})[["genome", "accession", "serovar"]]

# 4. Merge everything.
df = meta.merge(checkm, on="genome", how="left").merge(prokka, on="genome", how="left")
df["qc_verdict"] = df.apply(
    lambda r: "high quality"
    if pd.notna(r["Completeness"]) and pd.notna(r["Contamination"])
    and r["Completeness"] >= 90 and r["Contamination"] <= 5
    else "check",
    axis=1,
)
cols = ["genome", "accession", "serovar", "Completeness", "Contamination",
        "Genome_Size", "GC_Content", "cds", "trna", "rrna", "tmrna", "qc_verdict"]
cols = [c for c in cols if c in df.columns]
summary = df[cols].sort_values("Completeness", ascending=False)
summary.to_csv(os.path.join(RESULTS, "genome_qc_summary.csv"), index=False)
print(summary.to_string(index=False))
print("\nGenomes passing QC: %d / %d" % ((summary["qc_verdict"] == "high quality").sum(), len(summary)))

# 5. Figures.
plt.rcParams.update({"figure.dpi": 150, "font.size": 9})

# (a) Completeness vs contamination scatter, the classic QC plot.
fig, ax = plt.subplots(figsize=(7, 5))
colors = {"high quality": "#2ca02c", "check": "#d62728"}
for verdict, grp in summary.groupby("qc_verdict"):
    ax.scatter(grp["Contamination"], grp["Completeness"], s=80,
               c=colors.get(verdict, "grey"), label=verdict, edgecolors="black", linewidths=0.5)
for i, (_, r) in enumerate(summary.iterrows()):
    # vertical labels: all points share ~100% completeness, so horizontal
    # labels would collide in the tight contamination cluster
    ax.annotate(r["genome"], (r["Contamination"], r["Completeness"]),
                fontsize=7, rotation=90, ha="center", va="bottom",
                xytext=(0, 8), textcoords="offset points")
ax.axhline(90, color="grey", linestyle="--", linewidth=0.8)
ax.axvline(5, color="grey", linestyle="--", linewidth=0.8)
# zoom to the data: all genomes are near-perfect, the full 0-100 range
# would hide the real (small) differences between them
ax.set_xlim(-0.03, max(0.6, summary["Contamination"].max() * 1.3))
ax.set_ylim(99.7, 100.2)
ax.set_xlabel("CheckM2 contamination (%)")
ax.set_ylabel("CheckM2 completeness (%)")
ax.set_title("Genome completeness vs contamination (n=%d Salmonella assemblies)" % len(summary))
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(FIGURES, "01_completeness_contamination.png"))
plt.close(fig)

# (b) Annotated features per genome.
feat = summary.set_index("genome")[["cds", "trna", "rrna"]].dropna()
fig, ax = plt.subplots(figsize=(8, 4.5))
feat.plot(kind="bar", ax=ax, color=["#1f77b4", "#ff7f0e", "#2ca02c"], edgecolor="black", linewidth=0.3)
ax.set_ylabel("Count")
ax.set_title("Prokka-annotated features per genome")
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right", fontsize=8)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(FIGURES, "02_features_per_genome.png"))
plt.close(fig)

# (c) Genome size and GC content per genome.
fig, ax1 = plt.subplots(figsize=(8, 4.5))
x = range(len(summary))
ax1.bar(x, summary["Genome_Size"] / 1e6, color="#9467bd", edgecolor="black", linewidth=0.3)
ax1.set_ylabel("Genome size (Mb)", color="#9467bd")
ax1.set_xticks(list(x))
ax1.set_xticklabels(summary["genome"], rotation=30, ha="right", fontsize=8)
ax2 = ax1.twinx()
ax2.plot(x, summary["GC_Content"] * 100, color="#d62728", marker="o", linewidth=1.5, label="GC %")
ax2.set_ylabel("GC content (%)", color="#d62728")
ax1.set_title("Genome size and GC content per genome")
fig.tight_layout()
fig.savefig(os.path.join(FIGURES, "03_size_gc_per_genome.png"))
plt.close(fig)

print("Figures written to figures/.")
