import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

merged_path = snakemake.input[0]

fig_rats = snakemake.output["rat_sightings_per_zip"]
fig_avg_score = snakemake.output["avg_score_by_zip"]
fig_rat_crit = snakemake.output["rat_vs_critical"]
fig_corr_heatmap = snakemake.output["corr_heatmap"]
fig_rat_avg = snakemake.output["rat_vs_avg_score"]

merged_df = pd.read_csv(merged_path)

# 1. Rat sightings per ZIP
plt.figure(figsize=(10, 5))
df_sorted = merged_df.sort_values("rat_sightings", ascending=False)
plt.bar(df_sorted["Incident Zip"].astype(str), df_sorted["rat_sightings"])
plt.xticks(rotation=90)
plt.title("Rat Sightings per ZIP Code")
plt.xlabel("ZIP Code")
plt.ylabel("Number of Rat Sightings")
plt.tight_layout()
plt.savefig(fig_rats, bbox_inches="tight")
plt.close()

# 2. Average inspection score by ZIP
plt.figure(figsize=(10, 5))
plt.scatter(merged_df["Incident Zip"].astype(str), merged_df["avg_score"])
plt.xticks(rotation=90)
plt.title("Average Restaurant Inspection Score by ZIP Code")
plt.xlabel("ZIP Code")
plt.ylabel("Average Inspection Score")
plt.tight_layout()
plt.savefig(fig_avg_score, bbox_inches="tight")
plt.close()

# 3. Rat sightings vs critical violation rate
plt.figure(figsize=(6, 5))
sns.regplot(data=merged_df, x="critical_rate", y="rat_sightings")
plt.title("Rat Sightings vs Critical Violation Rate")
plt.tight_layout()
plt.savefig(fig_rat_crit, bbox_inches="tight")
plt.close()


# 4. Correlation heatmap
corr_vars = merged_df[
    ["rat_sightings", "avg_score", "median_score",
     "critical_rate", "total_inspections", "unique_restaurants"]
]
corr_matrix = corr_vars.corr()

plt.figure(figsize=(7, 5))
sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", fmt=".2f")
plt.title("Correlation Heatmap (ZIP-Level Aggregates)")
plt.tight_layout()
plt.savefig(fig_corr_heatmap, bbox_inches="tight")
plt.close()

# 5. Rat sightings vs average score
plt.figure(figsize=(6, 5))
sns.regplot(data=merged_df, x="avg_score", y="rat_sightings", scatter_kws={"alpha": 0.6})
plt.xlabel("Average Inspection Score")
plt.ylabel("Rat Sightings")
plt.title("Rat Sightings vs. Average Inspection Score")
plt.tight_layout()
plt.savefig(fig_rat_avg, bbox_inches="tight")
plt.close()
