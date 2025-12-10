# Snakefile

rule all:
    input:
        # cleaned restaurant and rat sightings datasets
        "data/clean/cleaned_rat.csv",
        "data/clean/cleaned_restaurant.csv",
        # merged + summary stats
        "data/merged/merged_df.csv",
        "results/summary_stats.csv",
        "results/corr_matrix.csv",
        # visual plots
        "figures/rat_sightings_per_zip.png",
        "figures/avg_score_by_zip.png",
        "figures/rat_vs_critical.png",
        "figures/corr_heatmap.png",
        "figures/rat_vs_avg_score.png"

# NOTE: the input "cleaned" datasets are exported from OpenRefine.
#       (they are NOT created by Snakemake).

rule unzip_rat:
    input:
        "data/clean/Rat-cleaned.zip"
    output:
        "data/clean/Rat-cleaned_in.csv"
    shell:
        "unzip -p {input} 'Rat-cleaned(in).csv' > {output}"

rule unzip_restaurant:
    input:
        "data/clean/NYC-cleaned.zip"
    output:
        "data/clean/NYC-cleaned_in.csv"
    shell:
        "unzip -p {input} 'NYC-cleaned(in).csv' > {output}"

rule select_clean_rat:
    input:
        "data/clean/Rat-cleaned_in.csv"
    output:
        "data/clean/cleaned_rat.csv"
    script:
        "scripts/select_clean_rat.py"

rule select_clean_restaurant:
    input:
        "data/clean/NYC-cleaned_in.csv"
    output:
        "data/clean/cleaned_restaurant.csv"
    script:
        "scripts/select_clean_restaurant.py"

rule merge_zip:
    input:
        rat="data/clean/cleaned_rat.csv",
        rest="data/clean/cleaned_restaurant.csv"
    output:
        "data/merged/merged_df.csv"
    script:
        "scripts/merge_zip.py"

rule summary_and_findings:
    input:
        "data/merged/merged_df.csv"
    output:
        summary="results/summary_stats.csv",
        corr="results/corr_matrix.csv"
    script:
        "scripts/summary_and_findings.py"

rule make_plots:
    input:
        "data/merged/merged_df.csv"
    output:
        rat_sightings_per_zip="figures/rat_sightings_per_zip.png",
        avg_score_by_zip="figures/avg_score_by_zip.png",
        rat_vs_critical="figures/rat_vs_critical.png",
        corr_heatmap="figures/corr_heatmap.png",
        rat_vs_avg_score="figures/rat_vs_avg_score.png"
    script:
        "scripts/make_plots.py"
