<h1>Final Project Report</h1>

<h2>Rat Sightings and Restaurant Inspection Outcomes in New York City</h2>

<h3>Contributors:</h3>
<ul>
  <li>Bushra Lat - blat2</li>
  <li>Sanvi Singh - sanvis2</li>
</ul>


<h3>Summary:</h3>


<p>Our project explores whether there is a meaningful relationship between rat sightings in New York City and the outcomes of restaurant health inspections across the same neighborhoods. Both issues are widely discussed in NYC, since rats are extremely common throughout the city and restaurant inspections are used to monitor public health. Since both datasets come from the same city and share the same geographic structure, we wanted to see whether neighborhoods that struggle more with rat sightings also tend to have restaurants with lower grades or more critical health violations. The goal of this project is not to claim that rats directly cause poor restaurant inspections, but instead to see whether these patterns appear together and what that might mean for public health and city planning.</p>
<p>We started by acquiring two open datasets from Kaggle. One dataset contains tens of thousands of reports of rat sightings made by the public. Each report includes the location of the sighting, including the ZIP code, the date the sighting was reported, and information about whether the complaint was inspected by the city. The second dataset contains restaurant inspection data from the NYC Department of Health. It includes the name, address, cuisine type, violation descriptions, inspection scores, and inspection grades for restaurants across all five boroughs. Because both datasets came from NYC Open Data originally, they follow similar organizational standards and are relatively compatible for integration after cleaning.</p>
<p>The general structure of our project followed the data lifecycle steps we learned in class. First, we acquired the data and documented the sources, licenses, and FAIR principles. Then we profiled and cleaned each dataset using OpenRefine. This step ended up being more intensive than expected. Many columns were inconsistent across both datasets, especially the address fields, business names, date formats, and ZIP codes. After cleaning, we selected the specific columns most relevant to our research questions and prepared the files for integration.</p>
<p>To join the two datasets, we used ZIP codes as the linking feature. Since both datasets include ZIP codes for each record, this made it possible to combine them in a meaningful way. However, because the restaurant dataset is much larger than the rat dataset, and because the number of restaurants is not evenly distributed across ZIP codes, we grouped each dataset by ZIP code before integrating them. This allowed us to create a final dataset where each row represents a single ZIP code and contains information about both rat sightings and restaurant inspection outcomes in that area.</p>
<p>Our analysis focuses on three main questions. First, we ask whether ZIP codes with higher numbers of rat sightings also have worse average inspection grades. Second, we ask whether higher rat-sighting density is connected to more critical violations in restaurants, especially ones related to rodents or food contamination. Third, we ask whether rodent-related violations specifically tend to be more common in areas with more rat activity. These questions help us understand whether the two datasets reflect similar patterns in the city’s neighborhoods.</p>
<p>By integrating the datasets and preparing them for analysis, we can now visualize rat activity alongside restaurant inspection outcomes. Although we have not completed our full statistical analysis yet, our early observations suggest that some ZIP codes with very high rat sightings also show higher inspection scores, meaning worse restaurant performance. This pattern is not universal, but it appears often enough that it seems worth investigating further. With the integrated dataset ready, we can proceed with measuring these relationships more formally and producing visualizations that can highlight trends more clearly.</p>
<p>Overall, this project demonstrates the value of integrating separate public datasets to explore broader public health questions. It also shows how much time and attention is required for cleaning and preparing real-world data before analysis can even begin. Our work so far sets up the next steps of the project, which include running exploratory visualizations, automating the workflow, and writing up a full reproducible report.</p>


<h3>Data Profile</h3>


<p>Our project uses two datasets that both come from public city information originally hosted on NYC Open Data. Although we downloaded them through Kaggle, both datasets were created by the NYC Department of Health and Mental Hygiene and the City of New York. Because they come from official city sources, the data is governed by open-data licensing, allowing for reuse and analysis as long as the original attribution is maintained. This makes them suitable for class projects and supports the FAIR principles of openness and reusability.</p>
<p>The first dataset contains records of rat sightings reported by the public. Each row represents one rat-related complaint submitted through the city’s 311 system. The dataset includes the date and time the complaint was created, the type of complaint, the location, the borough, the ZIP code, and the inspection status of the report. The dataset also includes latitude and longitude information, which makes it possible to map trends, although our project focuses on ZIP code because it is the variable that can connect the dataset to restaurant inspections. The dataset contains tens of thousands of observations covering many years. However, the volume of data also introduces challenges. Some records contain missing or incorrect ZIP codes, inconsistent date formats, and variations in how location fields were entered. We corrected these issues through data cleaning in OpenRefine.</p>
<p>The second dataset contains information about restaurant inspections conducted by the NYC Department of Health. Each row represents a single inspection event for a specific restaurant. The dataset includes restaurant names, addresses, cuisine types, inspection dates, violation codes, violation descriptions, inspection scores, and inspection grades. The scores represent how many points were assigned during the inspection, and the grades (A, B, or C) are determined based on the score. Restaurants with low scores receive better grades, while higher scores indicate more violations and result in worse grades. The dataset spans several years and includes hundreds of thousands of observations, making it much larger and more detailed than the rat dataset.</p>
<p>The restaurant dataset also has its own challenges. Because a restaurant can be inspected multiple times, the dataset contains many repeated business names and addresses. Business names appear in different formats, with abbreviations, misspellings, and inconsistent punctuation. Addresses have similar issues. These inconsistencies can affect grouping and analysis, so they require careful cleaning. In addition, some restaurants are missing grades or scores, which requires us to remove incomplete rows or handle missing data appropriately.</p>
<p>Both datasets contain information that is sensitive in the sense that it deals with public health, but the data does not include any personal identifiers. Rat complaints belong to properties, not individuals, and the restaurant dataset is already public information. Because of this, there are no major ethical concerns about the datasets. The information is intended to be used by the public and by analysts studying health and neighborhood conditions in NYC. The only important ethical responsibility for our project was to ensure that we handled the data according to the licensing rules and cited the sources correctly, which we have done.</p>
<p>To evaluate the datasets under the FAIR principles, we considered whether they are findable, accessible, interoperable, and reusable. Both datasets are findable through NYC Open Data and Kaggle and include metadata that explains their contents. They are accessible to anyone without requiring special permission or credentials. They are also interoperable because they contain clear variable names, simple formats, and common data types that make them easy to load in Python or OpenRefine. Finally, they are reusable because they include documentation, open licensing, and stable identifiers such as DOIs on Kaggle. Together, the datasets meet the FAIR principles and provide a strong foundation for integration.</p>
<p>Overall, the rat sightings dataset and restaurant inspections dataset each contain rich information about different aspects of NYC public health. When combined, they create a new opportunity to explore how environmental conditions and restaurant hygiene may relate to one another at the neighborhood level. The datasets have limitations, especially in terms of missing or inconsistent values, but after cleaning, they provide a reliable basis for our analysis.</p>


<h3>Data Quality</h3>


<p>Data quality played a major role in this project because both datasets contained many inconsistencies that needed to be resolved before integration. The rat sightings dataset and the restaurant inspections dataset were created by different city systems, which means they follow different patterns in how they record addresses, dates, names, and ZIP codes. We needed to examine the quality of the data carefully, identify problems, and fix them in a consistent way. Most of this work was done in OpenRefine, since it allowed us to explore patterns, cluster similar text values, and standardize fields efficiently.</p>
<p>One of the biggest issues in both datasets involved the address fields. The rat sightings dataset includes location fields such as street names, house numbers, and city names, but these are not always written consistently. For example, street names sometimes appear in full, such as “West 18th Street,” while other times they use abbreviations like “W 18 St.” Some entries contain repeated words or incorrectly formatted ordinal numbers. These issues make it hard to group observations by location or compare records. In the restaurant inspections dataset, the address issues were even more noticeable because restaurants appear multiple times across different inspections. Some restaurants had three or four versions of their address due to spelling differences or missing values. Cleaning these fields required clustering similar values, normalizing spacing and capitalization, and rewriting entries so that they follow one consistent format.</p>
<p>Another major area of concern was the restaurant names. Many restaurant names appear in multiple formats, especially chain restaurants. Dunkin Donuts and Baskin Robbins appeared in several different spellings or combinations. Local restaurants had even more inconsistent naming. To address this, we used OpenRefine’s clustering tools to group similar names together. When clusters did not provide a clear answer, we searched the restaurant online to confirm the correct spelling and selected that version as the standardized value. This step ensured that our grouped-by-ZIP-code aggregation would not overcount restaurants or produce misleading results.</p>
<p>ZIP codes were also a significant quality issue. Both datasets include ZIP codes, but many entries contained invalid or missing values. Some records had ZIP codes that did not exist, which could cause errors during integration. To fix this, we used the “us valid zips” list provided in class and used OpenRefine’s cross function to match each ZIP code to the official list. Records with invalid ZIP codes were removed, since they could not be reliably grouped or merged later. After cleaning, we confirmed that both datasets only contained valid, five-digit ZIP codes.</p>
<p>Date formats appeared in many different styles across both datasets. In the rat sightings dataset, dates included variations such as “04/12/2017,” “2017-04-12,” or even text-based forms. We standardized all dates into one consistent format using OpenRefine transformations. This correction made it possible to sort, filter, and analyze observations more easily. The restaurant dataset had fewer date problems, but some entries still required standardization.</p>
<p>Missing values were another issue we had to consider. While the datasets had thousands of complete entries, some columns had so many missing values that they were not useful for analysis. To handle this, we selected only the columns that were relevant to our research questions. After reducing each dataset to the necessary columns, most missing values disappeared. Any remaining missing entries were dropped using Python’s dropna function. This allowed us to create a clean and complete dataset for integration.</p>
<p>One final data quality concern was the imbalance between the two datasets. The restaurant dataset has several hundred thousand rows, while the rat dataset had a little over eighty thousand after cleaning. Because of this difference, integrating the datasets by row would not make sense. To address this, we aggregated each dataset by ZIP code before merging them. This created a more balanced structure and prevented the restaurant dataset from dominating the combined file.</p>
<p>Overall, the data cleaning and quality assessment steps were essential because they allowed us to create a reliable and consistent dataset for analysis. The corrections we made ensure that the results will be meaningful and that our workflow can be reproduced by others.</p>


<h3>Findings</h3>


<p>After cleaning, aggregating, and integrating the two datasets, we created a combined dataset where each row represents a single ZIP code. This allowed us to compare rat sightings and restaurant inspection outcomes in a clear and structured way. At this stage, our findings focus on the patterns we observed through our analysis and exploratory summaries. First, we looked at the summary statistics of all relevant variables in our integrated dataset. This included the variables for number of rat sightings per zip code, total number of restaurant inspections per zip code, number of unique restaurants per zip code, average restaurant inspection score per zip code, median restaurant inspection score per zip code, and average number of restaurants that were flagged as in “critical” condition per zip code, which we dubbed “critical rate”. The full Data Dictionary can be accessed <a href="https://github.com/bushralat/TheFAIRies/blob/main/Data%20Dictionary.md">here.</p>
<p>After looking at the summary statistics for each of these variables we were able to better understand how the data is distributed. For example, the mean for number of rat sightings per zip code was 450.4, but had a standard deviation of 453.96. This indicates a large variance in rat sightings per zip code. This is reflected in the figure <a href="https://github.com/bushralat/TheFAIRies/blob/main/rat_sightings_per_zip.png">“rat_sightings_per_zip”</a>, where you can see that there is a large variation in rat sightings per zip code. Certain clusters of zip codes, such as zip codes in the 11250 zip area, have thousands of rat sightings, whereas other zip codes have far fewer. This indicates that large amounts of rat sightings tend to be clustered in the same area and could be a sign of busier areas or areas with geographical issues.</p>
<p>In addition to looking at that, we also visualized the average restaurant inspection score by zip code, seen in the figure <a href = "https://github.com/bushralat/TheFAIRies/blob/main/avg_score_by_zip.png"> “avg_score_by_zip”</a>. This graph has individual plots for each average quality score per zip code. There is not a lot of variation in scores across all zip codes, but certain clusters of zip codes have more variation in average scores. There are also a few outliers with very good average scores (lower scores are better). This indicates that overall, restaurants score fairly consistently across zip codes with a few outliers.</p>
<p>After this, we looked at how rat sightings and restaurant inspections might be related. We created two scatter plots along with a line that represented a linear regression model fit. Each point in the scatter plot represents where the x and y-axis lay on the graph, while the line represents the overall trend of the graph. The first figure, <a href = "https://github.com/bushralat/TheFAIRies/blob/main/rat_vs_critical.png">“rat_vs_critical”</a>, depicts how the total number of right sightings and the critical rate of restaurants compare, per zip code. Within this graph, the points were scattered all over the graph, and did not follow a slope or the regression line. This indicates a weak relationship between rat sightings and how often restaurants are deemed “Critical”.</p>
<p>In addition to this, the figure <a href="https://github.com/bushralat/TheFAIRies/blob/main/rat_vs_avg_score.png">“rat_vs_avg_score”</a> depicts how rat sightings compare to the average restaurant score, per zip code. This plot had somewhat more of a trend between points on the graph, however the trend was still weak and did not indicate a strong correlation. Many of the zip codes had similarly high quality scores, but varied a lot in the number of rat sightings.</p>
<p>Finally, we visualized the numerical correlation between each variable in the figure <a href="https://github.com/bushralat/TheFAIRies/blob/main/corr_heatmap.png">“corr_heatmap”</a>. This figure states the correlation between each variable in our merged dataset, so we were able to see that there was not a high correlation between any of the restaurant inspection variables and the number of rat sightings by zip code. This confirmed our initial analysis that there is not a strong correlation between rat sightings and restaurant quality, at least as described by the datasets.</p>
<p>Overall, these visualizations and analysis tell us that high numbers of rat sightings in certain areas are not a strong indicator of restaurant quality. They could play a part in overall quality of living in certain areas of NYC, but do not have a strong relationship specifically with restaurant health and safety.</p>


<h3>Future Work</h3>


<p>Based on the analyses we made, we found that the answers to our initial questions came out to be inconclusive. We found no significant evidence that there is a strong relationship between rat sightings and restaurant quality in the area. When visualizing average restaurant quality score by zip code and number of rat sightings, we did see a weak positive correlation. This could be due to the fact that there is a weak relationship between rat sightings and restaurant scoring. However, it could also be because broader geographical issues are related to poor restaurant health and safety. For example, poorer areas in NYC could have more rat sightings and restaurants with more health code violations due to poor environmental conditions. In order to better understand the relationship between geographical conditions and restaurant quality, more data would need to be collected and analyzed. This data could be related to the environmental and geographic conditions of different neighborhoods in NYC. The data should include zip code so we can aggregate the datasets on zip code again. Then, we could explore areas such as neighborhood wealth, environmental conditions, and other aspects of neighborhood quality. We could then incorporate more granular geographical data analysis in order to explore relationships between geography and restaurant health quality.</p>
<p>Another aspect of our workflow that could be improved is within data cleaning. We did initial data cleaning of both datasets using OpenRefine, and then moved our data to a Jupyter Notebook to reduce the data to what we needed. We reduced the two datasets by dropping all variables except the ones we were interested in looking at. Furthermore, we dropped rows with missing values. Doing this dropped a chunk of our data. It wasn’t a majority of the datasets but it did reduce the data by a sizable amount. Therefore, we could approach data differently by implementing a pipeline that imputes new values in replacement of missing data rather than dropping all rows with missing values. This could improve our workflow by allowing us to utilize more of the data for our analysis of the integrated data.</p>
<p>Another thing we could implement into our workflow is further exploring the “Violations” variable within the restaurant inspections dataset. We dropped this variable from the restaurants inspections dataset before integrating the two datasets. However, if we kept this variable and categorized it based on the appearance of rodents in the violations, we could dive deeper into how rat sightings are related to rodent violations in NYC restaurants.</p>
<p>Finally, we could add more analysis to our workflow by creating a predictive model based on our data. This could be a linear regression model that predicts restaurant quality score based on rat sightings. Producing a predictive model based on our data could provide a little more insight as to how strong the relationship between rat sightings and restaurant quality is.</p>
<p>Overall, the ideas produced from this project could be further studied by acquiring new geographical data to explore new questions, and by making improvements to our workflow to better assess our analyses.</p>


<h3>Reproducing</h3>


Link to our Box File with raw and output data: https://uofi.box.com/s/f8dw62s6x4l1mobp91dv0fycrxc0hvku<a href="https://uofi.box.com/s/f8dw62s6x4l1mobp91dv0fycrxc0hvku">


<p>In order to reproduce our workflow, you will need the data after we cleaned it in OpenRefine.</p>
<p>This data is already in the <a href="https://github.com/bushralat/TheFAIRies/tree/main/data/clean">data/clean folder</a>.</p>
<p>The data inside the workflow is zipped due being large files, however the snakemake workflow unzips it.</p>

<p>The folder structure for the workflow is listed <a href="https://github.com/bushralat/TheFAIRies/blob/main/Snakemake%20Folder%20Structure">here</a>.</p>

<p>The following packages and virtual environment were used for this project:</p>
<ul>
  <li>Miniconda: ~/miniconda3/etc/profile.d/conda.sh</li>
  <li>Snakemake: conda install -c bioconda -c conda-forge snakemake</li>
  <li>Imports:</li>
      <li>Import pandas as pd</li>
      <li>Import matplotlib.pyplot as plt</li>
      <li>Import seaborn as sns</li>
</ul>


The workflow:

1. Starts from **OpenRefine–cleaned data** packaged in `.zip` files
2. Unzips and standardizes those cleaned files
3. Selects relevant columns and drops incomplete rows
4. Aggregates data by ZIP code and merges the two datasets (merged data found <a href="https://github.com/bushralat/TheFAIRies/tree/main/data/merged">here</a>.
5. Produces summary statistics and a correlation matrix (found in this <a href="https://github.com/bushralat/TheFAIRies/tree/main/results">folder</a>
6. Generates visual figures for analysis (found in this <a href="https://github.com/bushralat/TheFAIRies/tree/main/figures">folder</a>

## Script Descriptions

### `select_clean_rat.py`
Cleans the OpenRefine cleaned rat sightings dataset by selecting relevant columns, removing incomplete rows, printing basic data quality checks, and saving a fully cleaned CSV (cleaned_rat.csv).

### `select_clean_restaurant.py`
Cleans the OpenRefine-processed restaurant inspection dataset by selecting key fields, dropping rows with missing data, and producing a cleaned dataset (cleaned_restaurant.csv).

### `merge_zip.py`
Aggregates the cleaned rat and restaurant datasets by ZIP code, computes statistics such as rat sighting counts, inspection totals, and critical violation rates, merges the results, and outputs a unified ZIP-level dataset. Produces "merged_df".

### `summary_and_findings.py` *(or `summary_and_corr.py`)*
Generates numerical summary statistics and a correlation matrix from the merged ZIP-level dataset and saves the resulting CSV files into the `results/` directory.

### `make_plots.py`
Creates all figures used in the analysis—including bar charts, scatterplots, regression plots, histograms, and a correlation heatmap—and saves them as PNG files in the `figures/` directory.

---

## Snakemake Unzip Rules

### `unzip_rat`
A Snakemake rule that extracts the rat sightings CSV from the zipped OpenRefine-cleaned file, creating the input needed for rat data cleaning.

### `unzip_restaurant`
A Snakemake rule that extracts the restaurant inspection CSV from the zipped OpenRefine-cleaned file, creating the input needed for restaurant data cleaning.

### `Snakefile`
Defines the full automated workflow, including unzipping, cleaning, merging, summarizing, and plotting. Snakemake determines appropriate execution order and ensures full reproducibility.


<h3>References</h3>


<p>Original Sources:</p>
<p>New York City Department of Health and Mental Hygiene. (n.d.). Restaurant Inspection Results [Data set]. NYC Open Data. https://data.cityofnewyork.us/Health/Restaurant-Inspection-Results/43nn-pn8j</p>
<p>New York City 311 / NYC Department of Health and Mental Hygiene. (n.d.). Rodent Inspection [Data set]. NYC Open Data. https://data.cityofnewyork.us/Health/Rodent-Inspection/p937-wjvj</p>

<p>Redistributed Versions Used for Download:</p>
<p>Boysen, J. (2017). NYC Restaurant Inspections [Data set]. Kaggle. https://www.kaggle.com/datasets/new-york-city/nyc-inspections</p>
<p>Boysen, J. (2017). NYC Rat Sightings [Data set]. Kaggle. https://www.kaggle.com/datasets/new-york-city/nyc-rat-sightings/data</p>
<p>OpenRefine Team. (2023). OpenRefine [Software]. https://openrefine.org</p>
<p>Microsoft. (2023). Visual Studio Code [Software]. https://code.visualstudio.com</p>


<h3>Metadata Documentation</h3>


<p><b>DataCite-Style MetaData for Project Datasets</b></p>
<p><b>Primary Dataset 1: NYC Restaurant Inspection Results</b></p>

<p><u>Identifier:</u></p>
<p>Official Source URL: https://data.cityofnewyork.us/Health/Restaurant-Inspection-Results/43nn-pn8j</p>

<p><u>Title:</u></p>
<p>Restaurant Inspection Results</p>

<p><u>Creator / Publisher:</u></p>
<p>New York City Department of Health and Mental Hygiene (DOHMH)</p>

<p><u>Publication Year:</u></p>
<p>Ongoing dataset, accessed 2025</p>

<p><u>Resource Type:</u></p>
<p>Open government dataset (CSV)</p>

<p><u>Description:</u></p>
<p>This dataset contains detailed records of restaurant inspections for food establishments across New York City. It includes inspection dates, inspection types, violation descriptions, violation codes, risk categories, cuisines, restaurant names, and the health grades assigned by DOHMH inspectors. The records span multiple years and provide insight into restaurant sanitation, food safety compliance, and health code violations at the ZIP-code level. These variables allow us to analyze location-based inspection trends and compare them with rat sighting patterns.</p>

<p><u>Version / Access Details:</u></p>
<p>Accessed through NYC Open Data in 2025. Downloaded through Kaggle, which mirrors the official NYC dataset but the authoritative provenance is NYC Open Data.</p>

<p><u>License:</u></p>
<p>Open Data Commons Public Domain Dedication and License (PDDL) — standard for NYC Open Data.</p>

<p><u>Related Identifiers:</u></p>
<p>Kaggle mirror used for download: https://www.kaggle.com/datasets/new-york-city/nyc-inspections</p>


<p><b>Primary Dataset 2: NYC Rodent Inspection (Rat Sightings) Dataset</b></p>

<p><u>Identifier:</u></p>
<p>Official Source URL: https://data.cityofnewyork.us/Health/Rodent-Inspection/p937-wjvj</p>

<p><u>Title:</u></p>
<p>Rodent Inspection (Rat Sightings)</p>

<p><u>Creator/Publisher:</u></p>
<p>New York City Department of Health and Mental Hygiene (DOHMH)</p>

<p><u>Publication Year:</u></p>
<p>Ongoing dataset, accessed 2025</p>

<p><u>Resource Type:</u></p>
<p>Open government dataset (CSV)</p>

<p><u>Description:</u></p>
<p>This dataset contains records of rodent-related complaints, inspections, and outcomes reported across New York City. Each record includes the date of the sighting or inspection, the complaint type, inspection result (Passed, Failed, Monitoring, etc.), the geographic coordinates, and the ZIP code. These records allow us to analyze the distribution of rat activity in different areas of NYC and compare those patterns to restaurant inspection performance.</p>

<p><u>Version / Access Details:</u></p>
<p>Accessed through NYC Open Data in 2025. Downloaded through Kaggle, which redistributed an older version of the same data. The authoritative source is NYC Open Data.</p>

<p><u>License:</u></p>
<p>Open Data Commons Public Domain Dedication and License (PDDL).</p>

<p><u>Related Identifiers:</u></p>
<p>Kaggle mirror used for download: https://www.kaggle.com/datasets/new-york-city/nyc-rat-sightings/data</p>


<p><b>Integrated Dataset (Created for This Project)</b></p>

<p><u>Identifier:</u></p>
<p>Created by: Sanvi Singh and Bushra Lat, Not publicly released outside this repository</p>

<p><u>Title:</u></p>
<p>Integrated Dataset Combining NYC Rat Sightings and Restaurant Inspection Results</p>

<p><u>Creator/Publisher:</u></p>
<p>Sanvi Singh and Bushra Lat</p>

<p><u>Publication Year:</u></p>
<p>2025</p>

<p><u>Resource Type:</u></p>
<p>Derived Dataset</p>

<p><u>Description:</u></p>
<p>This dataset was created by merging the cleaned NYC Rat Sightings dataset and the NYC Restaurant Inspection Results dataset using ZIP codes as the shared key. Before integration, both datasets were cleaned using OpenRefine to standardize address formats, correct business names, validate ZIP codes, and normalize date formats. Missing or invalid ZIP codes were removed, and irrelevant variables were filtered out. The merged dataset includes rat sighting counts, restaurant inspection counts, aggregated inspection scores, the number of unique restaurants in each ZIP code, and the percentage of restaurants earning an “A” grade. This dataset supports analysis of the relationship between rat activity and restaurant sanitation outcomes across NYC.</p>

<p><u>License:</u></p>
<p>Derived from NYC public domain datasets. This dataset inherits NYC Open Data’s PDDL license.</p>


<h3>Software and Tools Used</h3>


<p><b>OpenRefine:</b></p>
<p>Version: 3.x</p>
<p>Description: Used to clean and standardize both datasets, including address normalization, ZIP code validation, removing duplicates, merging restaurant name variants, converting case formats, and harmonizing date formats.</p>
<p>URL: https://openrefine.org<a href="https://openrefine.org"></p>


<p><b>Visual Studio Code</b></p>
<p>Version: 1.106.3</p>
<p>Description: Used to run Python scripts for grouping, aggregating, merging, and exporting the integrated dataset.</p>
<p>URL: https://code.visualstudio.com<a href="https://code.visualstudio.com"></p>


<p><b>Python + Pandas</b></p>
<p>Used for grouping by ZIP code, aggregating fields, dropping missing values, and merging datasets.</p>
