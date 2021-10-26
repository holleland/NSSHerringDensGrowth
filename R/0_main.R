# -------- Main script --------- 
# -- Author: Sondre Hølleland --
# ------------------------------
# The following scripts will fill up the folders:
# - data
# - tables
# - plots
# which is the output used in the paper. 

# Set up the data
source("R/1_data.R")
# Make figures only based on the data: 
source("R/2_data_figures.R")
# Estimate VBGF year-by-year
source("R/3_estimate_year_by_year.R")
# Once this is done, we can run all the models:
source("R/4_model_runs.R")
# - Produces figure S2-S3, S6-S7
# Produce figures: 
source("R/5_result_figures.R")
# Produce tables: 
source("R/6_result_tables.R")
# Length vs weight figures (and model results)
source("R/7_length_to_weight.R")
