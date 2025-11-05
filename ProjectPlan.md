# IS 477 Course Project Plan

## Overview
This project will explore the relationship between rat sightings in NYC and the details of restaurant inspections, focusing on possible links to rat sightings and the safety/quality of restaurants in NYC. We will use two open datasets from Kaggle: one on restaurant inspections in NYC and another on rat sightings in NYC. Our main goal is to find patterns that show how rat sightings may correlate with poor grades from restaurant inspections.

The project will combine both datasets to see if zip codes with high levels of rat sightings, especially ones that are not taken care of, can be connected to zip codes with restaurants that fail health inspections. We hope that by analyzing these two sources together, we can better understand how restaurants can maintain healthier conditions and rat sightings can be reduced in NYC. The project will also follow the full data lifecycle from collecting and organizing data to cleaning, integrating, and analyzing it. Each step will be documented carefully to make the work reproducible and transparent.

## Research Questions
Our work will focus on three main questions:

- Do ZIP codes with higher numbers of rat sightings also report worse restaurant inspection grades?
- Is there a correlation between rat sighting density and the frequency of critical health violations in restaurants (e.g., vermin, food storage, contamination)?
- Are rodent related violations in restaurant inspections more common in areas with more rat sightings?

## Team
- **Sanvi Singh** – Responsible for collecting the data, checking ethical requirements, documenting metadata, and doing early data exploration.  
- **Bushra Lat** – Responsible for cleaning and integrating the data, automating the workflow, and creating visualizations.  

## Datasets

### NYC Rat Sightings Dataset
**Source:** [Kaggle NYC Rat Sightings](https://www.kaggle.com/datasets/new-york-city/nyc-rat-sightings/data)

This dataset contains information about rat sightings in NYC, with each row being an individual rat sighting observed. As rats are very prominent in NYC, and were observed to be about ¼ of the human population in NYC in 2014, people can file a complaint when they see a rat in private and public properties. Then, the rat incident is expected to be handled by the party responsible. This dataset contains details about location of each sighting and the date it occurred, as well as how each incident was handled among thousands of observations.

### NYC Restaurant Inspections
**Source:** [Kaggle NYC Restaurant Inspections](https://www.kaggle.com/datasets/new-york-city/nyc-inspections)

This dataset contains information about restaurant inspections for permitted food establishments in NYC. Restaurants are graded on an A-F scale during inspections done by the city’s health department. The dataset includes descriptors such as address, cuisine description, inspection date, type, action, violation code and description(s). The data covers all of NYC and starts Jan 1, 2010 up until Aug 29, 2017. The data was collected by the NYC Department of Health.

## Integration Plan
The datasets will be integrated using zip code as the common variable. Before merging, we will standardize variable names and values so that the data can be compared properly. Integration will be done using Python and pandas. We will document all steps clearly to make the process reproducible.

## Timeline
1. **Week 1:** We will download and explore both datasets and document their sources and licenses.  
2. **Week 2:** We will review the structure and quality of the data, noting missing or inconsistent values.  
3. **Week 3:** We will clean and organize the data so that it can be analyzed easily.  
4. **Week 4:** We will integrate them and begin running exploratory analyses and visualizations to identify any connections.  
5. **Final Weeks:** We will automate the workflow, write up our findings, and prepare the final report for submission.  

## Constraints
One limitation of this project is that the datasets come from different studies, so certain zip codes may be covered more than others. Additionally, the NYC Restaurant Inspections data is a much bigger dataset and has many more observations compared to the NYC Rat Sightings dataset. This means that our data may be unbalanced and preprocessing would be required to balance the two datasets so they can be joined. Another challenge is that the datasets have different formats and naming conventions, which will require some cleaning and standardization before integration. Finally, we will need to handle missing values through preprocessing and confirm that both datasets comply with open data terms.

## Gaps and Next Steps
Some gaps within our datasets include disparities in dataset size. The Restaurant Inspections dataset is much bigger than the Rat Sightings dataset, so we will need to take a smaller section of data points from the Restaurant Inspection’s dataset to integrate with the Rat Sightings dataset. This could affect our analyses of the combined datasets, and preprocessing will be required in order to integrate the two datasets.

Another gap which requires further research is deciding which variables from both datasets are most relevant to each other. Both datasets include features about zip code/location for each rat sighting and restaurant inspection. We would also like to determine which features are most indicative of rat infestations and poor restaurant health safety when the two datasets are integrated. The Rat Sightings dataset alone has 50+ variables so we will have to determine which variables from that dataset are most relevant to our questions. We will address this gap by doing research on the individual variables within both datasets, and doing preprocessing and exploratory analysis of both datasets in order to determine which variables are most relevant to answering our research questions.

## Connection to Course Concepts
Our project is connected to what we are learning in the course about the data lifecycle. From acquiring the dataset to assessing it for ethics, we are applying what we learned about the FAIR data principles and data acquisition so that we can analyze the data using Python. Additionally, we are applying what we have learned in data modeling in order to assess our selected feature variables and how they affect our predictor variables. This will allow us to address our research questions about correlations between restaurant inspections and rat sightings. Furthermore, we will apply what we have learned in data integration in order to integrate our two datasets and cross-analyze the different features in both datasets. This will allow us to understand the relationship between restaurant health safety and rat prevalence.

## References
- Jacob Boysen. (2017). NYC Restaurant Inspections [Data set]. Kaggle. 
  https://www.kaggle.com/datasets/new-york-city/nyc-inspections

- Jacob Boysen. (2017). NYC Rat Sightings [Data set]. Kaggle.   
  https://www.kaggle.com/datasets/new-york-city/nyc-rat-sightings/data
