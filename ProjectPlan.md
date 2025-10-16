# IS 477 Course Project Plan

## Overview
This project will explore the relationship between sleep disorders and cognitive decline, focusing on possible links to Alzheimer’s disease. We will use two open datasets from Kaggle: one on sleep disorder diagnoses and another on Alzheimer’s disease. Our main goal is to find patterns that show how sleep quality, physical health, and daily habits may relate to changes in cognitive health.

The project will combine both datasets to see if lifestyle factors like sleep duration or stress levels can be connected to signs of dementia or memory loss. We hope that by analyzing these two sources together, we can better understand how sleep and mental health affect long-term brain health. The project will also follow the full data lifecycle from collecting and organizing data to cleaning, integrating, and analyzing it. Each step will be documented carefully to make the work reproducible and transparent.

## Research Questions
Our work will focus on three main questions:

- Is there a relationship between sleep disorder factors such as sleep duration or stress level and the cognitive performance indicators used to measure Alzheimer’s disease?
- Can sleep-related variables help predict early signs of Alzheimer’s disease or cognitive decline?
- How do age, gender, and health factors influence the connection between sleep and brain health?

## Team
- **Sanvi Singh** – Responsible for collecting the data, checking ethical requirements, documenting metadata, and doing early data exploration.  
- **Bushra Lat** – Responsible for cleaning and integrating the data, automating the workflow, and creating visualizations.  

## Datasets

### Sleep Disorder Diagnosis Dataset
**Source:** [Kaggle Sleep Disorder Diagnosis Dataset](https://www.kaggle.com/datasets/mdsultanulislamovi/sleep-disorder-diagnosis-dataset)

This dataset includes information about individuals’ sleep patterns and physical health. It contains details such as sleep duration, quality, stress level, heart rate, body mass index, occupation, and whether a sleep disorder like insomnia or sleep apnea is present. It is a CSV file with around 400 rows. We will use variables such as age, gender, sleep duration, and disorder type to understand how they relate to overall health and behavior. The dataset is anonymized and publicly available for educational use.

### Alzheimer’s Disease Dataset
**Source:** [Kaggle Alzheimer’s Disease Dataset](https://www.kaggle.com/datasets/rabieelkharoua/alzheimers-disease-dataset)

This dataset contains information about patients with different stages of Alzheimer’s disease, ranging from nondemented to fully demented. It includes demographic and clinical variables such as age, gender, education level, memory test scores, and MRI-based brain volume measures. It is also a CSV file. We will use this dataset to study patterns in cognitive decline and link them with lifestyle or physical health factors from the sleep dataset. The dataset is open for research and contains no personal identifiers.

## Integration Plan
The datasets will be integrated using age and gender as the common variables. Before merging, we will standardize variable names and values so that the data can be compared properly. Integration will be done using Python and pandas. We will document all steps clearly to make the process reproducible.

## Timeline
1. **Week 1:** We will download and explore both datasets and document their sources and licenses.  
2. **Week 2:** We will review the structure and quality of the data, noting missing or inconsistent values.  
3. **Week 3:** We will clean and organize the data so that it can be analyzed easily.  
4. **Week 4:** We will integrate them and begin running exploratory analyses and visualizations to identify any connections.  
5. **Final Weeks:** We will automate the workflow, write up our findings, and prepare the final report for submission.  

## Constraints
One limitation of this project is that the datasets come from different studies and do not contain the same individuals. This means that our findings will show associations rather than direct cause-and-effect relationships. Another challenge is that the datasets have different formats and naming conventions, which will require some cleaning and standardization before integration. The sample sizes are also fairly small, which might limit the strength of the results. Finally, we will need to handle missing values and confirm that both datasets comply with open data terms.  

## Gaps and Next Steps
Some gaps within our datasets include disparities in dataset size. The Alzheimer’s dataset is much bigger than the Sleep Disorder Diagnosis dataset, so we will need to take a smaller section of data points from the Alzheimer’s dataset to integrate with the Sleep Disorder Diagnosis dataset. This could affect our analyses of the combined datasets, and preprocessing will be required in order to integrate the two datasets.

Another gap which requires further research is deciding which variables from both datasets are most relevant to each other. Both datasets include features about health history for each individual patient, and our goal is to determine which features are most indicative of sleep disorders and Alzheimer’s respective to each dataset and when the two datasets are integrated. We will address this gap by doing research on the individual variables within both datasets, and doing preprocessing and exploratory analysis of both datasets in order to determine which variables are most relevant to answering our research questions.

## Connection to Course Concepts
Our project is connected to what we are learning in the course about the data lifecycle. From acquiring the dataset to assessing it for ethics, we are applying what we learned about the FAIR data principles and data acquisition so that we can analyze the data using Python. Additionally, we are applying what we have learned in data modeling in order to assess our selected feature variables and how they affect our predictor variables. This will allow us to address our research questions about what facts impact sleep and cognitive health. Furthermore, we will apply what we have learned in data integration in order to integrate our two datasets and cross-analyze the different features in both datasets. This will allow us to understand the relationship between sleep health and cognitive health.

## References
- Sultanul, O. (2025). *Sleep Disorder Diagnosis Dataset* [Data set]. Kaggle.  
  https://www.kaggle.com/datasets/mdsultanulislamovi/sleep-disorder-diagnosis-dataset  

- El Kharoua, R. (2024). *Alzheimer’s Disease Dataset* [Data set]. Kaggle.  
  https://www.kaggle.com/datasets/rabieelkharoua/alzheimers-disease-dataset
