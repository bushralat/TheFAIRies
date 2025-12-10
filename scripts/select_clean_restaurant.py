import pandas as pd

rest_input = snakemake.input[0]    # "data/clean/NYC-cleaned(in).csv"
rest_output = snakemake.output[0]  # "data/clean/cleaned_restaurant.csv"

cleaned_restaurant_df = pd.read_csv(rest_input)

# check missing values
print("Restaurant (OpenRefine cleaned) null counts:")
print(cleaned_restaurant_df.isnull().sum())

cleaned_restaurant_df = cleaned_restaurant_df[
    [
        'CAMIS','DBA','BORO','BUILDING','STREET','ZIPCODE','PHONE',
        'CUISINE DESCRIPTION','INSPECTION DATE','ACTION','VIOLATION CODE',
        'VIOLATION DESCRIPTION','CRITICAL FLAG','SCORE','RECORD DATE',
        'INSPECTION TYPE'
    ]
]

# drop missing values
df_restaurant = cleaned_restaurant_df.dropna()

# dtypes of each column
print("Restaurant dtypes after column selection and dropna:")
print(df_restaurant.dtypes)

# final cleaned restaurant inspections file
df_restaurant.to_csv(rest_output, index=False)
