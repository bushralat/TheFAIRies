import pandas as pd

# input for cleaned restaurant and rat sightings files
rat_path = snakemake.input["rat"]   # data/clean/cleaned_rat.csv
rest_path = snakemake.input["rest"] # data/clean/cleaned_restaurant.csv
merged_output = snakemake.output[0]

# create final cleaned dataframes
df_rat = pd.read_csv(rat_path)
df_restaurant = pd.read_csv(rest_path)

# group restaurant dataset by zipcode for merging
restaurant_grouped = (
    df_restaurant.groupby("ZIPCODE")
    .agg(
        total_inspections=("CAMIS", "count"),
        unique_restaurants=("CAMIS", "nunique"),
        avg_score=("SCORE", "mean"),
        median_score=("SCORE", "median"),
        critical_rate=("CRITICAL FLAG", lambda x: (x == "Critical").mean())
    )
    .reset_index()
)

# group rat dataset by zipcode for merging
rat_grouped = (
    df_rat
        .groupby("Incident Zip")
        .agg(
            rat_sightings=("Unique Key", "count")
        )
        .reset_index()
)

# merge restaurant and rat sightings dataset on zip
merged_df = rat_grouped.merge(
    restaurant_grouped,
    left_on="Incident Zip",
    right_on="ZIPCODE",
    how="left"
)

# fill in missing values within merged dataset with 0
merged_df = merged_df.fillna(0)

# create CSV file of merged dataframe
merged_df.to_csv(merged_output, index=False)