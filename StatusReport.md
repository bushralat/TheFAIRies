<h2>Project Overview</h2>

<p>For this project, we are working on integrating a dataset detailing observations of NYC restaurant inspections and a dataset detailing observations of NYC rat sightings. The purpose of this project is built around the following research questions:</p>

<ul>
  <li>Do ZIP codes with higher numbers of rat sightings also report worse restaurant inspection grades?</li>
  <li>Is there a correlation between rat sighting density and the frequency of critical health violations in restaurants (e.g., vermin, food storage, contamination)?</li>
  <li>Are rodent related violations in restaurant inspections more common in areas with more rat sightings?</li>
</ul>

<p>The tasks for this project include:</p>
<ul>
  <li>Downloading and exploring both datasets, as well as documenting their sources.</li>
  <li>Reviewing the structure and quality of the data, noting missing or inconsistent values.</li>
  <li>Cleaning and organizing the data so that it can be analyzed easily.</li>
  <li>Integrating them and running exploratory analyses and visualizations to identify any connections.</li>
  <li>Automating the workflow.</li>
</ul>

<h2>Current Progress</h2>

<p>So far for the project, we have acquired the datasets that we will use. We acquired the data from Kaggle, which was sourced from the NYC Department of Health and the City of New York, from NYC Open Data. We downloaded these datasets directly from Kaggle. We evaluated the data based on the FAIR guiding principles to determine if the data is findable, accessible, interoperable, and reusable. We determined that both datasets score well within the FAIR guiding principles because the data is publically available from the city of New York, has DOIs, is free to download, uses common vocabulary for integration, and contains information about metadata. We used OpenRefine to clean the datasets. After cleaning the data using OpenRefine, we picked the specific columns we wanted to look at within each dataset and integrated the two datasets using Python.</p>

<p>The supporting artifacts for our progress include:</p>
<ul>
  <li><a href="https://github.com/bushralat/TheFAIRies/blob/main/open_refine_history_rat_sightings.json">OpenRefine history for the rat sightings file.</a></li>
  <li><a href="https://github.com/bushralat/TheFAIRies/blob/main/open_refine_history_NYC_Restaurant.json">OpenRefine history for the restaurant inspections file.</a></li>
  <li><a href="https://github.com/bushralat/TheFAIRies/blob/main/cleaned_rat.csv.zip">The cleaned rat sightings dataset.</a></li>
  <li><a href="https://github.com/bushralat/TheFAIRies/blob/main/cleaned_restaurant.csv.zip">The cleaned restaurant sightings dataset.</a></li>
  <li><a href="https://github.com/bushralat/TheFAIRies/blob/main/merged_df.csv">The integrated dataset.</a></li>
</ul>

<h2>Updated Timeline</h2>

<p>Within the next two weeks we will run exploratory analyses on the integrated dataset and identify possible connections. We will also produce an automated end-to-end workflow and produce metadata and data documentation to support understandability and reusability.</p>

<h2>Changes to Project Plan</h2>

<p>We made a major change to our project after Milestone 2. Our original plan was to use two Kaggle datasets related to sleep disorders and Alzheimer’s disease, and we planned to integrate them using age and gender. However, after starting the early exploration stage, we realized that the two health datasets did not match well for integration. They came from completely different studies, had very different structures, and did not contain compatible variables beyond basic demographics. Because of this, it became clear that merging them would not produce meaningful or reliable results. The datasets were also much smaller than we expected, and the sleep dataset in particular did not have enough depth to support our research questions about cognitive decline. After receiving feedback suggesting we choose two datasets with a stronger natural connection, we decided to switch topics.</p>

<p>Our new plan focuses on NYC rat sightings and NYC restaurant inspections. These datasets are both from NYC Open Data and are much more suitable for integration because they share a common geographic variable: ZIP code. This makes it possible to analyze them together in a meaningful way. The new topic also allowed us to form clearer research questions about how rat activity in neighborhoods relates to restaurant health violations. Another change we made is that the new datasets are much larger and required more detailed cleaning, especially the address fields. This shifted part of our workload and timeline, since more time had to be spent on data cleaning and standardizing the two files before integration. Even though the topic changed, our workflow still follows the data lifecycle we originally planned, including acquisition, cleaning, integration, analysis, and automation.</p>

<p>The main change to our project plan was replacing the original Sleep and Alzheimer’s datasets with NYC rat sightings and NYC restaurant inspection data. This change made the project more feasible and gave us a clearer path to answering our research questions. The rest of our plan stayed mostly the same, but we updated our methods and timeline to match the new datasets and the amount of cleaning they required.</p>

<h2>Team Member Contributions</h2>

<h3>Sanvi’s contributions:</h3>

<p>For this milestone, I spent most of my time cleaning and standardizing both datasets, and this ended up being a lot more work than I expected. The hardest part was fixing the street and address fields because the same locations were written in so many different ways. Streets had different abbreviations, random spacing, misspellings, and repeated words. For example, “West 18 Street” might also show up as “W 18 Street,” “West 18 Street,” or even “600 West West 161 Street.” I fixed these using OpenRefine by normalizing directions, correcting spacing, converting everything to title case, and rewriting ordinal numbers like “Second” or “3” so they were consistent. I also cleaned the business names in the restaurant dataset, which was another big task. Many restaurants appeared under several different spellings, and I had to merge them into one consistent version. Some of these were confusing enough that I had to look them up online to make sure I was choosing the correct spelling. For example, Dunkin Donuts and Baskin Robbins appeared in many different formats, and some smaller restaurants were spelled in a way that didn’t seem right, so I did a bit of research just to confirm the proper name before standardizing it.</p>

<p>Another major part of my work involved validating ZIP codes in both datasets. I used the “us valid zips” file we were given and used a cell.cross() function in OpenRefine to check each ZIP code against the official list. After matching them, I removed all rows with invalid or missing ZIP codes so the datasets would join correctly later. I also cleaned up the date and time formats in the rat sightings dataset, since they were written in different styles. I converted everything into one standard date format to make sorting and filtering easier. On top of all that, I standardized many text columns by converting them to title case so the data would be cleaner and easier to read.</p>

<p>Overall, the cleaning process took a long time because there were so many small inconsistencies across thousands of rows. Even though it was time-consuming, it makes the datasets much more reliable for the next steps of the project. I feel like the cleaning work I completed puts us in a strong position to begin merging and analyzing the data.<br>
To explore the datasets, we exported the data as CSV files and used VSCode to analyze the data with Python.</p>

<h3>Bushra’s Contributions:</h3>

<p>After the datasets were cleaned in OpenRefine, I uploaded the datasets into VSCode so we could integrate the two datasets. First, I checked if there were still missing values in each of the datasets. A few of the columns in both datasets had a considerable amount of missing values. I picked the columns that were relevant to our research questions and created new datasets that only included these columns. This got rid of a large amount of missing values. A couple of the columns still had missing values, so I used .dropna to drop all the remaining missing values. This got rid of all missing values left. After that I integrated the two datasets.</p>

<p>In order to integrate the two datasets, I grouped each dataset by zip code. I did this because the restaurant inspection dataset had 373,848 rows after cleaning, and the rat dataset had 82,424 rows after cleaning. I grouped each dataset by zip code so that each row in the integrated dataset was a unique zip code. I chose specific columns to aggregate into the dataset. In addition to incident ZIP, those columns include: Incident Zip, rat_sightings, ZIPCODE, total_inspections, unique_restaurants, avg_score, median_score, and grade_A_pct. These columns allow us to look at rat sightings as well as restaurant quality in specific locations (defined by zip code). I also manually filled in missing zip codes with 0, because some zip codes may exist in one dataset but not the other. This resulted in our integrated dataset, titled “merged_df”.</p>
