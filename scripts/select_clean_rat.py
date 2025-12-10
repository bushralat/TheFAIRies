import pandas as pd

rat_input = snakemake.input[0]   # "data/clean/Rat-cleaned(in).csv"
rat_output = snakemake.output[0] # "data/clean/cleaned_rat.csv"

cleaned_rat_df = pd.read_csv(rat_input)

# check missing values
print("Rat (OpenRefine cleaned) null counts:")
print(cleaned_rat_df.isnull().sum())

cleaned_rat_df = cleaned_rat_df[
    [
        'Unique Key','Created Date','Closed Date','Agency','Agency Name',
        'Complaint Type','Descriptor','Location Type','Incident Zip',
        'Incident Address','Street Name','Address Type','City','Status',
        'Due Date','Resolution Action Updated Date','Community Board',
        'Borough','Location'
    ]
]

# drop missing values
df_rat = cleaned_rat_df.dropna()

# dtypes of each column
print("Rat dtypes after column selection and dropna:")
print(df_rat.dtypes)

# final cleaned rat sightings file
df_rat.to_csv(rat_output, index=False)
