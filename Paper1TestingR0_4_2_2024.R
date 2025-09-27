####Paper1TestingR0_4_2_2024 ####
#### Load necessary packages####
library(tidyverse)
library(reshape2)
##load regression packages##
library(tidyverse)  # data wrangling and visualization
library(knitr)      # beautifying tables
library(car)        # for checking assumptions, e.g. vif etc.
library(broom)      # for tidy model output
library(questionr)  # for odds.ratios
library(sjPlot)     # for plotting results of log.regr.
library(sjmisc)     # for plotting results of log.regr.
library(effects)    # for probability output and plots
library(emmeans)
require(sandwich)
require(msm)
library(dplyr)
library(purrr)
library(tibble)


##load graphics packages##
library(latticeExtra) 
library(RColorBrewer)
library(plotly)

#load EpiEstim Package
library(EpiEstim)

# Set working directory (replace with the path to your folder)
setwd("/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024")
####Load Data#####
load("Analyticdf_paper1_freq0.rdata")

##We can run estimate_R on the incidence data to estimate the reproduction number R. 
#For this, we need to specify i) the time window(s) over which to estimate R and 
#ii) information on the distribution of the serial interval.
#For i), the default behavior is to estimate R over weekly sliding windows. 
#This can be changed through the config$t_start and config$t_end arguments 
#(see below, “Changing the time windows for estimation”). 
#For ii), there are several options, specified in the method argument.
#The simplest is the parametric_si method, where you only specify the mean and standard deviation of the SI.

#Define a function that computes estimate_R for each simulation in week 1:
compute_estimate_R_wk1 <- function(data){
  require(EpiEstim)  #the function estimate_R is from the EpiEstim package
  
  result <- estimate_R(incid = data$Incidence, 
                       method = "parametric_si", 
                       config = make_config(list(mean_si = 5.2, std_si = 0.15)))
  return((result$R[1, , drop = FALSE]))
}

# This will give you a dataframe where each row corresponds to a unique simulation number, 
# and there is a list-column called estimate_results which contains the results from estimate_R for each simulation.
results_list_wk1 <- analyticdf %>%
  group_by(unique_sim_number) %>%
  nest() %>%
  mutate(estimate_results = map(data, compute_estimate_R_wk1))

results_df_wk1 <- results_list_wk1 %>%
  unnest(cols = c(estimate_results))

meanR0wk1<- mean(results_df_wk1$`Mean(R)`) #mean R0 is 1.81
medianR0wk1<-median(results_df_wk1$`Median(R)`) #median R0 is 1.39


#Define a function that computes estimate_R for each simulation in week 2:
compute_estimate_R_wk2 <- function(data){
  require(EpiEstim)  #the function estimate_R is from the EpiEstim package
  
  result <- estimate_R(incid = data$Incidence, 
                       method = "parametric_si", 
                       config = make_config(list(mean_si = 5.2, std_si = 0.15)))
  return((result$R[2, , drop = FALSE]))
}

# This will give you a dataframe where each row corresponds to a unique simulation number, 
# and there is a list-column called estimate_results which contains the results from estimate_R for each simulation.
results_list_wk2 <- analyticdf %>%
  group_by(unique_sim_number) %>%
  nest() %>%
  mutate(estimate_results = map(data, compute_estimate_R_wk2))

results_df_wk2 <- results_list_wk2 %>%
  unnest(cols = c(estimate_results))

meanR0wk2<-mean(results_df_wk2$`Mean(R)`) 
medianR0wk2<-median(results_df_wk2$`Median(R)`) 

R0table_tau0.0125<- cbind(meanR0wk1,medianR0wk1, meanR0wk2, medianR0wk2)
save(R0table_tau0.0125, file ="R0table_tau0.0125.Rdata")

####Estimate R for all weeks####

compute_estimate_R <- function(data){
  require(EpiEstim)
  
  result <- estimate_R(incid = data$Incidence, 
                       method = "parametric_si", 
                       config = make_config(list(mean_si = 5.2, std_si = 0.15)))
  
  # Assuming 'R' is a matrix or dataframe within the result object
  # and you're interested in the first 15 weeks
  return(result$R[1:15, , drop = FALSE])
}

results_list <- analyticdf %>%
  group_by(unique_sim_number) %>%
  nest() %>%
  mutate(estimate_results = map(data, compute_estimate_R))

# Create an empty dataframe to store the final results
final_results_df <- data.frame()

# Loop through each row in results_list
for(i in seq_len(nrow(results_list))) {
  
  sim_number <- results_list$unique_sim_number[i]
  
  # Extract the 15 rows of R values for the i-th simulation
  r_values_15_weeks <- results_list$estimate_results[[i]]
  
  # Add a new column for unique_sim_number
  r_values_15_weeks$unique_sim_number <- sim_number
  
  # Append to the final dataframe
  final_results_df <- rbind(final_results_df, r_values_15_weeks)
}

# Now, final_results_df will contain the R values for each of the 15 weeks,
# and each set of 15 rows will have the same unique_sim_number.
