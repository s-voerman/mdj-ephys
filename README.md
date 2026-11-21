# mdj-ephys
Repository containing analysis scripts for the Journal of Physiology paper 'Electrophysiological properties of mesodiencephalic junction neurons projecting to the inferior olive'
A link to the article can be found here:
https://physoc.onlinelibrary.wiley.com/doi/10.1113/JP291269 

All scripts are written in Python 3.14+

# Initial statement
This is the initial and minimal version of the scripts that were used for data analysis and figure generation. This repository is currently still a work-in-progress. I currently  don't have much time to maintain/add to this repo. I will do my best to supply everything that is reasonable to include ASAP. Nothing in this repository has been tested to actually run, I will try to get to a version that will run on processed data (if raw data is not easily available) so that all scripts can be verified more easily. I've tried to prioritize the scripts that are important to verify, 

My apologies for the current incompleteness of this projects scripts. If you have specific requests or questions about any data or scripts in the paper feel free to raise an issue on this page or contact the corresponding author. 

This is what is currently included:
1. Analysis scripts in `analysis`. This folder contains the primary scripts used in data analysis
2. Processed data in `data`
3. Cell metadata in `metadata`

# Data & metadata
The `data` folder contains processed data. I'm currently looking into options of also providing raw data, but the file sizes may be too large to allow for that. As an in between I've supplied the outputs of several scripts in `analysis`. Im trying to figure out currently how to include IPSC data, as this is a very large dataset that can't really be reduced in size as much as the other data.

The `metadata` folder contains details about the cells which were analyzed

# Analysis
The `analysis` folder contains analysis scripts. These transform the raw data into the data found in the `data` folder. Subfolders contain:
1. `ephys-main` which contains most ephys analysis scripts and a script for the PCA and k-means clustering
2. `ipsc` which contains scripts used in automated analysis of IPSCs and PTX IPSC experiments
3. `optogenetics` which contains scripts used in optogenetic analysis

# Figures
The `figures` folder has a similar strcuture to the `analysis` folder and is also split in the same categories. 

# Statistics
The `statistics` folder contains scripts that were used in the statistical analysis of data and primarily uses processed data in the `data` folder

# Output
The `output` folder will be for outputting plot pngs/pdfs and data analysis outputs

# Issues & questions
For questions please contact the corresponding author of the article.