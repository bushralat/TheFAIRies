import pandas as pd

merged_path = snakemake.input[0]
summary_output = snakemake.output["summary"]
corr_output = snakemake.output["corr"]

merged_df = pd.read_csv(merged_path)

summary_df = merged_df[
    [
        'rat_sightings',
        'total_inspections',
        'unique_restaurants',
        'avg_score',
        'median_score',
        'critical_rate'
    ]
].describe()

summary_df.to_csv(summary_output)

corr_features = merged_df[
    ["rat_sightings", "avg_score", "median_score",
     "critical_rate", "total_inspections", "unique_restaurants"]
]
corr_matrix = corr_features.corr()
corr_matrix.to_csv(corr_output)
