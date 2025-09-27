####Paper_1_Analysis####

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

##load graphics packages##
library(latticeExtra) 
library(RColorBrewer)
library(plotly)
#load EpiEstim
library(EpiEstim)


# Set working directory (replace with the path to your folder)
setwd("/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024")

####Load manalyticdf Data#####

# List all files with the pattern 'Manalyticdf_paper1_'
files <- list.files(pattern = "Manalyticdf_paper1_freq\\d+\\.RData")

# Initialize an empty list to store the dataframes
dfs_list <- list()

# Loop through all files and load them
for (i in 1:length(files)) {
  file <- files[i]
  
  # Load the .RData file
  load(file)
  
  # Rename the dataframe
  assign(paste0("minianalyticdf_", i), minianalyticdf)
  
  # Add the renamed dataframe to the list
  dfs_list[[i]] <- get(paste0("minianalyticdf_", i))
  
  # Remove the original dataframe from the environment
  rm(minianalyticdf)
}
# Combine all dataframes in the list into a single dataframe
minianalyticdf <- bind_rows(dfs_list)

# Loop through the numbers 1 to 8 and remove the dataframes
for (i in 1:8) {
  rm(list = paste0("minianalyticdf_", i))
}
#Convert frequency scale from daily to weekly
minianalyticdf$frequency<- minianalyticdf$frequency*7
minianalyticdf$freq<- minianalyticdf$freq*7

########Load analyticdf Data#####
# List all files with the pattern 'Analyticdf_paper3_'
files <- list.files(pattern = "Analyticdf_paper1_freq\\d+\\.RData")

# Initialize an empty list to store the dataframes
dfs_list <- list()

# Loop through all files and load them
for (i in 1:length(files)) {
  file <- files[i]
  
  # Load the .RData file
  load(file)
  
  # Rename the dataframe
  assign(paste0("Analyticdf_", i), analyticdf)
  
  # Add the renamed dataframe to the list
  dfs_list[[i]] <- get(paste0("Analyticdf_", i))
  
  # Remove the original dataframe from the environment
  rm(analyticdf)
}
# Combine all dataframes in the list into a single dataframe
analyticdf <- bind_rows(dfs_list)

# Loop through the numbers 1 to 8 and remove the dataframes
for (i in 1:8) {
  rm(list = paste0("Analyticdf_", i))
}
#remove dfs_list
rm(dfs_list)
#Convert frequency scale from daily to weekly
analyticdf$frequency<- analyticdf$freq
analyticdf$frequency<- analyticdf$frequency*7
analyticdf$freq<- analyticdf$freq*7
# Convert Time to numeric if it isn't already
analyticdf$Time <- as.numeric(analyticdf$Time)

#calculate cumulative incidence
analyticdf <- analyticdf %>%
  arrange(unique_sim_number, Time, frequency) %>%
  group_by(unique_sim_number, frequency) %>%
  mutate(cumulative_incidence = cumsum(Incidence))

#calculate Outbreak20 variable
analyticdf$Outbreak20<- ifelse(analyticdf$everinfected>=20,"yes", "no")

####Test if still growing at 100 days####

# Function to fit exponential model if growth rate is positive
is_exponentially_growing <- function(data) {
  # Filter to the last 20 days of data for the simulation
  last_20_days <- data %>% filter(Time > 80)
  
  # Perform a log transformation of Incidence to fit a linear model
  last_20_days <- last_20_days %>%
    mutate(log_incidence = log(Incidence + 1)) # Adding 1 to avoid log(0)
  
  # Fit a linear model
  fit <- lm(log_incidence ~ Time, data = last_20_days)
  
  # Extract the coefficient of Time (which corresponds to b in the exponential model)
  growth_rate <- coef(fit)["Time"]
  
  # Return TRUE if the growth rate is positive
  return(growth_rate > 0)
}

# Apply the exponential growth check to each simulation where Outbreak50 = yes
df_growth_exp_50 <- analyticdf %>%
  filter(Outbreak50 == "yes") %>%
  group_by(unique_sim_number) %>% 
  do(is_growing = is_exponentially_growing(.)) %>%
  summarise(proportion_still_growing_exp = mean(is_growing))

#calculate proportion of simulations that are growing 
prop_growing_e_50<-sum(df_growth_exp_50)/nrow(df_growth_exp_50)

# Apply the exponential growth check to each simulation where Outbreak20 = yes
df_growth_exp_20 <- analyticdf %>%
  filter(Outbreak20 == "yes") %>%
  group_by(unique_sim_number) %>% 
  do(is_growing = is_exponentially_growing(.)) %>%
  summarise(proportion_still_growing_exp = mean(is_growing))

#calculate proportion of simulations that are growing 
prop_growing_e_20<-sum(df_growth_exp_20)/nrow(df_growth_exp_20)

# Apply the exponential growth check to each simulation where Outbreak10 = yes
df_growth_exp_10 <- analyticdf %>%
  filter(Outbreak10 == "yes") %>%
  group_by(unique_sim_number) %>% 
  do(is_growing = is_exponentially_growing(.)) %>%
  summarise(proportion_still_growing_exp = mean(is_growing))

#calculate proportion of simulations that are growing 
prop_growing_e_10<-sum(df_growth_exp_10)/nrow(df_growth_exp_10)

#combine into table
prop_growing_100days<-cbind(prop_growing_e_50, prop_growing_e_20, prop_growing_e_10)
prop_growing_100days<- as.data.frame(prop_growing_100days)
names(prop_growing_100days)<- c("Proportion of outbreaks of 50+ growing at 100 days", "Proportion of outbreaks of 20+ growing at 100 days", "Proportion of outbreaks of 10+ growing at 100 days")
write.csv(prop_growing_100days, file = "prop_outbreaks_still_growing.csv")


####Estimate Contributions to FOI####
#Load parameter table, except for parameters that are already in analyticdf
par_tab <- expand.grid(tauL=0.0125, epsilonL=0.19, 
                       delta=0.77, gammaL=0.1, piL=0.02, tauG=0.0125, 
                       epsilonG=0.19,thetaG=0, omegaG=0, 
                       chiG=0.01, gammaG=0.1,piG=0, epsilonIS=0.19, 
                       gammaIS=0.1, epsilonIAL = 0.0396, 
                       gammaIAL=0.07, epsilonIAG=0.0396, gammaIAG=0.07, 
                       iota=0.2, psi=0.75, zeta=0.5, mu=0.15, nu=0.48)
#add par_tab as new variables to analyticdf
analyticdf<-cbind(analyticdf,par_tab)
#Calculate effective number of infected neighbors per time step
analyticdf$effNinfneighbors<-analyticdf$ISL + analyticdf$IPL*analyticdf$psi +analyticdf$ISG +analyticdf$IPG*analyticdf$psi +analyticdf$IAL*analyticdf$zeta +analyticdf$IAG*analyticdf$zeta +analyticdf$IPISO*analyticdf$psi*analyticdf$iota + analyticdf$ISISO*analyticdf$iota

##calculate FOI:(probability of infection of any susceptible per time step)
#analyticdf$FOI<- 1-(1-analyticdf$tauL)^analyticdf$Ninfneighbors
##calculate share of FOI from each disease state for each time step
analyticdf$shareFOI_IPL<- ifelse(analyticdf$effNinfneighbors>0,analyticdf$IPL*analyticdf$psi / analyticdf$effNinfneighbors,0) 
analyticdf$shareFOI_IPG<- ifelse(analyticdf$effNinfneighbors>0,analyticdf$IPG*analyticdf$psi / analyticdf$effNinfneighbors,0) 
analyticdf$shareFOI_ISL<- ifelse(analyticdf$effNinfneighbors>0,analyticdf$ISL / analyticdf$effNinfneighbors,0) 
analyticdf$shareFOI_ISG<- ifelse(analyticdf$effNinfneighbors>0,analyticdf$ISG / analyticdf$effNinfneighbors,0) 
analyticdf$shareFOI_IAL<- ifelse(analyticdf$effNinfneighbors>0,analyticdf$IAL*analyticdf$zeta/ analyticdf$effNinfneighbors,0) 
analyticdf$shareFOI_IAG<- ifelse(analyticdf$effNinfneighbors>0,analyticdf$IAG*analyticdf$zeta / analyticdf$effNinfneighbors,0) 
analyticdf$shareFOI_IPISO<- ifelse(analyticdf$effNinfneighbors>0,analyticdf$IPISO*analyticdf$psi*analyticdf$iota / analyticdf$effNinfneighbors,0) 
analyticdf$shareFOI_ISISO<- ifelse(analyticdf$effNinfneighbors>0,analyticdf$ISISO*analyticdf$iota / analyticdf$effNinfneighbors,0) 
analyticdf$sumofFOIshares<- ifelse(analyticdf$effNinfneighbors>0,analyticdf$shareFOI_IPL+analyticdf$shareFOI_IPG+analyticdf$shareFOI_ISL+analyticdf$shareFOI_ISG+ analyticdf$shareFOI_IAL + analyticdf$shareFOI_IAG+analyticdf$shareFOI_IPISO+ analyticdf$shareFOI_ISISO,0) 
##calculate number of effective infected neighbors for each disease state for each time step
analyticdf$effNinfIPL<-analyticdf$shareFOI_IPL*analyticdf$effNinfneighbors
analyticdf$effNinfIPG<- analyticdf$shareFOI_IPG*analyticdf$effNinfneighbors
analyticdf$effNinfISL<- analyticdf$shareFOI_ISL*analyticdf$effNinfneighbors
analyticdf$effNinfISG<-analyticdf$shareFOI_ISG*analyticdf$effNinfneighbors
analyticdf$effNinfIAL<-analyticdf$shareFOI_IAL*analyticdf$effNinfneighbors
analyticdf$effNinfIAG<- analyticdf$shareFOI_IAG*analyticdf$effNinfneighbors
analyticdf$effNinfIPISO<-analyticdf$shareFOI_IPISO*analyticdf$effNinfneighbors
analyticdf$effNinfISISO<-analyticdf$shareFOI_ISISO*analyticdf$effNinfneighbors

#Calculate the total effective infectious time for each state for each simulation (across all time steps).
analyticdf <- analyticdf %>%
  group_by(unique_sim_number) %>%
  mutate(
    totalEffNinfIPL = sum(effNinfIPL, na.rm = TRUE),
    totalEffNinfIPG = sum(effNinfIPG, na.rm = TRUE),
    totalEffNinfISL = sum(effNinfISL, na.rm = TRUE),
    totalEffNinfISG = sum(effNinfISG, na.rm = TRUE),
    totalEffNinfIAL = sum(effNinfIAL, na.rm = TRUE),
    totalEffNinfIAG = sum(effNinfIAG, na.rm = TRUE),
    totalEffNinfIPISO = sum(effNinfIPISO, na.rm = TRUE),
    totalEffNinfISISO = sum(effNinfISISO, na.rm = TRUE),
    totalEffNinfNeighbors = sum(effNinfneighbors, na.rm = TRUE)
  ) %>%
  ungroup()
#Calculate the proportion of each state's total effective infectious
#time to the overall total effective infectious time across all states for 
#each simulation.

analyticdf <- analyticdf %>%
  mutate(
    propIPL = totalEffNinfIPL / totalEffNinfNeighbors,
    propIPG = totalEffNinfIPG / totalEffNinfNeighbors,
    propISL = totalEffNinfISL / totalEffNinfNeighbors,
    propISG = totalEffNinfISG / totalEffNinfNeighbors,
    propIAL = totalEffNinfIAL / totalEffNinfNeighbors,
    propIAG = totalEffNinfIAG / totalEffNinfNeighbors,
    propIPISO = totalEffNinfIPISO / totalEffNinfNeighbors,
    propISISO = totalEffNinfISISO / totalEffNinfNeighbors
  )

#get the average contribution of each state across all simulations, 
analyticdf <- analyticdf %>%
  mutate(
    avgPropIPL = mean(propIPL, na.rm = TRUE),
    avgPropIPG = mean(propIPG, na.rm = TRUE),
    avgPropISL = mean(propISL, na.rm = TRUE),
    avgPropISG = mean(propISG, na.rm = TRUE),
    avgPropIAL = mean(propIAL, na.rm = TRUE),
    avgPropIAG = mean(propIAG, na.rm = TRUE),
    avgPropIPISO = mean(propIPISO, na.rm = TRUE),
    avgPropISISO = mean(propISISO, na.rm = TRUE)
  )

###Compare the proportions of effective infectious time spent in each state 
# between simulations with different cumulative incidences at Time=100

# Extract rows where Time=100
time_100_df <- analyticdf %>% filter(Time == 100)

# Identify the unique simulation numbers for each condition
high_50_sim_numbers <- time_100_df %>% 
  filter(cumulative_incidence >= 50) %>% 
  pull(unique_sim_number)

low_49_sim_numbers <- time_100_df %>% 
  filter(cumulative_incidence <= 49) %>% 
  pull(unique_sim_number)

# Calculating average proportions for identified high_50_sim_numbers
high_50_avg_props <- analyticdf %>% 
  filter(unique_sim_number %in% high_50_sim_numbers) %>% 
  summarise(across(starts_with("prop"), mean, na.rm = TRUE))
high_50_avg_props<- round(high_50_avg_props,3)
# Calculating average proportions for identified low_49_sim_numbers
low_49_avg_props <- analyticdf %>% 
  filter(unique_sim_number %in% low_49_sim_numbers) %>% 
  summarise(across(starts_with("prop"), mean, na.rm = TRUE))
low_49_avg_props<- round(low_49_avg_props,3)

# Calculate the grand mean of prop variables across all simulations
grand_mean_props <- analyticdf %>%
  summarise(across(starts_with("prop"), mean, na.rm = TRUE))
grand_mean_props<- round(grand_mean_props,3)

# Create combined table of proportions for sims above and below 50: props_by_epidemic_size
high_50_avg_props$EpidemicSize <- "50 or More Cumulative Infections"
low_49_avg_props$EpidemicSize <- "Fewer than 50 Cumulative Infections"
grand_mean_props$EpidemicSize <- "All Epidemics"

props_by_epidemic_size<- rbind(high_50_avg_props,low_49_avg_props, grand_mean_props)

props_by_epidemic_size<- props_by_epidemic_size[ ,c("EpidemicSize","propIPL","propIPG","propISL","propISG","propIAL","propIAG","propIPISO","propISISO")]
names(props_by_epidemic_size)<- c("Epidemic Size","Infectious, Pre-Symptomatic Lab Workers","Infectious, Pre-Symptomatic General Population","Infectious, Symptomatic Lab Workers","Infectious, Symptomatic General Population","Infectious, Asymptomatic Lab Workers","Infectious, Asymptomatic General Population","Infectious, Pre-Symptomatic, Isolated","Infectious, Symptomatic, Isolated")
props_by_epidemic_size<- transpose(props_by_epidemic_size)
rownames(props_by_epidemic_size)<-c("Epidemic Size","Infectious, Pre-Symptomatic Lab Workers","Infectious, Pre-Symptomatic General Population","Infectious, Symptomatic Lab Workers","Infectious, Symptomatic General Population","Infectious, Asymptomatic Lab Workers","Infectious, Asymptomatic General Population","Infectious, Pre-Symptomatic, Isolated","Infectious, Symptomatic, Isolated")
# save eff_inf_time_props_by_epidemic_size.csv
write.csv(props_by_epidemic_size,"eff_inf_time_props_by_epidemic_size.csv", row.names = TRUE)



# Calculate the mean of prop variables grouped by frequency
frequency_mean_props <- analyticdf %>%
  group_by(frequency) %>%
  summarise(across(starts_with("prop"), mean, na.rm = TRUE))
frequency_mean_props<-round(frequency_mean_props,3)

frequency_mean_df <- analyticdf %>%
  group_by(frequency) %>%
  summarise(
    AvgPropIPL = mean(propIPL, na.rm = TRUE),
    AvgPropIPG = mean(propIPG, na.rm = TRUE),
    AvgPropISL = mean(propISL, na.rm = TRUE),
    AvgPropISG = mean(propISG, na.rm = TRUE),
    AvgPropIAL = mean(propIAL, na.rm = TRUE),
    AvgPropIAG = mean(propIAG, na.rm = TRUE),
    AvgPropIPISO = mean(propIPISO, na.rm = TRUE),
    AvgPropISISO = mean(propISISO, na.rm = TRUE)
  )
frequency_mean_df<-round(frequency_mean_df,3)

#



####**Create Epidemic Plots**####

####Plot Inc Sample of Sims Daily####

# Step 1: Filter the DataFrame for frequency values of 0 and 5
filtered_analyticdf <- analyticdf %>% 
  filter(frequency %in% c(0, 5))

# Step 2: Randomly sample 3 unique simulations for each frequency value
set.seed(42) # For reproducibility

# For frequency = 0
sampled_sims_0 <- filtered_analyticdf %>% 
  filter(frequency == 0) %>%
  select(unique_sim_number) %>%
  distinct() %>%
  sample_n(min(3, n())) %>%
  pull(unique_sim_number)

# For frequency = 5
sampled_sims_5 <- filtered_analyticdf %>% 
  filter(frequency == 5) %>%
  select(unique_sim_number) %>%
  distinct() %>%
  sample_n(min(3, n())) %>%
  pull(unique_sim_number)

# Combine the sampled simulation numbers
sampled_sims <- c(sampled_sims_0, sampled_sims_5)

# Step 3: Use the filtered dataframe in ggplot2
sampled_analyticdf <- filtered_analyticdf %>% 
  filter(unique_sim_number %in% sampled_sims)

# Note that sampled_analyticdf should now be much smaller and quicker to plot
ggplot(sampled_analyticdf, aes(x = Time, y = Incidence)) +
  geom_line(aes(color = as.factor(unique_sim_number), group = unique_sim_number)) +
  facet_wrap(~ frequency) +
  theme_minimal() +
  labs(title = "Incidence over Time for Sampled Simulations",
       x = "Time",
       y = "Incidence")

####Plot Inc sample of sims Weekly####
# Step 0: Create a week variable
analyticdf <- analyticdf %>% 
  mutate(Week = ceiling(Time / 7))

# Step 1: Filter the DataFrame for frequency values of 0 and 5
filtered_analyticdf <- analyticdf %>% 
  filter(frequency %in% c(0, 5))

# Step 2: Aggregate by week, simulation, and frequency
weekly_data <- filtered_analyticdf %>% 
  group_by(Week, unique_sim_number, frequency) %>% 
  summarise(Weekly_Incidence = sum(Incidence)) %>% 
  ungroup()

# Step 3: Randomly sample 5 unique simulations for each frequency value
set.seed(1) # For reproducibility

# For frequency = 0
sampled_sims_0 <- weekly_data %>% 
  filter(frequency == 0) %>%
  select(unique_sim_number) %>%
  distinct() %>%
  sample_n(min(5, n())) %>%
  pull(unique_sim_number)

# For frequency = 5
sampled_sims_5 <- weekly_data %>% 
  filter(frequency == 5) %>%
  select(unique_sim_number) %>%
  distinct() %>%
  sample_n(min(5, n())) %>%
  pull(unique_sim_number)

# Combine the sampled simulation numbers
sampled_sims <- c(sampled_sims_0, sampled_sims_5)

# Step 4: Use the filtered dataframe in ggplot2
sampled_weekly_data <- weekly_data %>% 
  filter(unique_sim_number %in% sampled_sims)

# Plot plot_sampled_weekly_incidence_by_freq
plot_sampled_weekly_incidence_by_freq<- ggplot(sampled_weekly_data, aes(x = Week, y = Weekly_Incidence)) +
  geom_line(aes(color = as.factor(unique_sim_number), group = unique_sim_number)) +
  facet_wrap(~ frequency, labeller = as_labeller(c(`0` = "0 Tests Per Week", `5` = "5 Tests Per Week"))) +
  labs(title = "Weekly Incidence Over Time For Sampled Simulations",
       x = "Week",
       y = "Weekly Incidence") +
  scale_color_discrete(name = "Unique Simulation Number") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold")) # This line adjusts the facet label text
plot_sampled_weekly_incidence_by_freq
ggsave("plot_sampled_weekly_incidence_by_freq.png", plot = plot_sampled_weekly_incidence_by_freq, width = 10, height = 6, dpi = 300)


####Plot Cum Inc sample of sims weekly####

# Compute cumulative incidence for each simulation
sampled_weekly_data <- sampled_weekly_data %>%
  group_by(unique_sim_number, frequency) %>%
  arrange(Week) %>%
  mutate(Cumulative_Incidence = cumsum(Weekly_Incidence)) %>%
  ungroup()
# Plot plot_sampled_weekly_cumulative_incidence_by_freq
plot_sampled_weekly_cumulative_incidence_by_freq<-ggplot(sampled_weekly_data, aes(x = Week, y = Cumulative_Incidence)) +
  geom_line(aes(color = as.factor(unique_sim_number), group = unique_sim_number)) +
  facet_wrap(~ frequency, labeller = as_labeller(c(`0` = "0 Tests Per Week", `5` = "5 Tests Per Week"))) +
  scale_color_discrete(name = "Unique Simulation Number") +
  labs(title = "Cumulative Incidence over Time for Sampled Simulations",
       x = "Week",
       y = "Cumulative Incidence") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold")) +
  scale_color_discrete(name = "Unique Simulation Number")

plot_sampled_weekly_cumulative_incidence_by_freq
ggsave("plot_sampled_weekly_cumulative_incidence_by_freq.png", plot = plot_sampled_weekly_cumulative_incidence_by_freq, width = 10, height = 6, dpi = 300)

####Plot Outbreak Sizes####
# Filter the dataframe for Time = 100 
data_at_100 <- analyticdf %>%
  filter(Time == 100)

# Calculate total number of simulations
total_sims <- nrow(data_at_100)

# Plot the histogram
plot_outbreak_sizes<-ggplot(data_at_100, aes(x = cumulative_incidence)) +
  geom_histogram(breaks = seq(0, max(data_at_100$cumulative_incidence, na.rm = TRUE) + 10, by = 5), 
                 fill = "blue", 
                 color = "black",
                 alpha = 0.7,
                 aes(y = ..count../total_sims * 100)) +  # Convert counts to percentages
  scale_x_continuous(breaks = seq(0, 300, by = 10)) +
  labs(title = "Distribution of Cumulative Incidence at 100 Days",
       x = "Cumulative Incidence at 100 Days",
       y = "Percentage of Simulations") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # To display y-axis values as percentages
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 8),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"))

plot_outbreak_sizes
ggsave("plot_outbreak_sizes.png", plot = plot_outbreak_sizes, width = 10, height = 6, dpi = 300)



####Plot Mean incidence by time####

# Convert Time to numeric if it isn't already
analyticdf$Time <- as.numeric(analyticdf$Time)

# Create a new variable for week classification
analyticdf$Week <- (analyticdf$Time - 1) %/% 7 + 1

#Calculate the Weekly Incidence for Each Simulation and Week:
analyticdf_weekly <- analyticdf %>%
  group_by(unique_sim_number, Week) %>%
  summarise(weekly_incidence = sum(Incidence, na.rm = TRUE)) %>%
  ungroup()
#Calculate the Mean Weekly Incidence for Each Week:
summary_weekly <- analyticdf_weekly %>%
  group_by(Week) %>%
  summarise(mean_weekly_incidence = mean(weekly_incidence, na.rm = TRUE),
            sd_weekly_incidence = sd(weekly_incidence, na.rm = TRUE),
            se_weekly_incidence = sd_weekly_incidence / sqrt(n()),
            upper = mean_weekly_incidence + se_weekly_incidence * 1.96,
            lower = mean_weekly_incidence - se_weekly_incidence * 1.96) %>%
  arrange(Week)

##Plot the Mean Weekly Incidence##

# Initial plot without X and Y axis
plot(summary_weekly$Week, summary_weekly$mean_weekly_incidence, type = "n", 
     xlab = "Week", ylab = "Mean Weekly Incidence",
     ylim = c(0,2),
     main = "Mean Weekly Incidence by Week", xaxt='n', yaxt='n')
# Add the lines for mean and CI
lines(summary_weekly$Week, summary_weekly$mean_weekly_incidence, col = "blue")
lines(summary_weekly$Week, summary_weekly$upper, col = "red", lty = 2)
lines(summary_weekly$Week, summary_weekly$lower, col = "red", lty = 2)
# Add custom tick marks for X-axis
axis(1, at = unique(summary_weekly$Week), labels = unique(summary_weekly$Week))
# Add custom tick marks for Y-axis between 0 and 2 in increments of 0.5
axis(2, at = seq(0, 2, by = 0.5), labels = seq(0, 2, by = 0.5))
# Add legend
legend("topright", legend = c("Mean", "95% CI"), col = c("blue", "red"), lty = c(1, 2))

####Mean weekly incidence by frequency####

#Calculate the Weekly Incidence for Each Simulation, Week, and Frequency value
analyticdf_weekly <- analyticdf %>%
  group_by(unique_sim_number, Week, frequency) %>%
  summarise(weekly_incidence = sum(Incidence, na.rm = TRUE))

# Calculate the mean and standard deviation of weekly incidence for each Week and Frequency value
mean_weekly_by_frequency <- analyticdf_weekly %>%
  group_by(Week, frequency) %>%
  summarise(
    mean_weekly_incidence = mean(weekly_incidence, na.rm = TRUE),
    sd_weekly_incidency = sd(weekly_incidence, na.rm = TRUE),
    se_weekly_incidence = sd(weekly_incidence, na.rm = TRUE) / sqrt(n())
  )

# Plot the mean weekly incidence by Frequency value
plot_mean_weekly_incidence_by_freq<-ggplot(mean_weekly_by_frequency, aes(x = Week, y = mean_weekly_incidence, color = as.factor(frequency))) +
  geom_line() +
  scale_x_continuous(breaks = seq(min(mean_weekly_by_frequency$Week), max(mean_weekly_by_frequency$Week), by = 1)) +
  labs(title = "Mean Weekly Incidence by Test Frequency Per Week",
       x = "Week",
       y = "Mean Weekly Incidence",
       color = 
         "Weekly
Test 
Frequency") +
  theme_minimal() +
  theme(
    axis.ticks.x = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold")
  )
plot_mean_weekly_incidence_by_freq
ggsave("plot_mean_weekly_incidence_by_freq.png", plot_mean_weekly_incidence_by_freq)

# # Plot the mean weekly incidence by Frequency value with standard deviation
# ggplot(mean_sd_weekly_by_frequency, aes(x = Week, y = mean_weekly_incidence, color = as.factor(frequency))) +
#   geom_line() +
#   geom_ribbon(aes(
#     ymin = mean_weekly_incidence - sd_weekly_incidency,
#     ymax = mean_weekly_incidence + sd_weekly_incidency,
#     fill = as.factor(frequency)), alpha = 0.2) +
#   scale_x_continuous(breaks = seq(min(mean_sd_weekly_by_frequency$Week), max(mean_sd_weekly_by_frequency$Week), by = 1)) +
#   labs(title = "Mean Weekly Incidence by Frequency with SE",
#        x = "Week",
#        y = "Mean Weekly Incidence",
#        color = "Weekly
# Test 
# Frequency",
#        fill = "frequency") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.ticks.x = element_line(color = "black"),
#     legend.position = "bottom"
#   )

####Mean Cumu Incidence by Frequency####
# Calculate cumulative incidence for each simulation, time, and frequency
analyticdf_cumulative <- analyticdf %>%
  arrange(unique_sim_number, Time, frequency) %>%
  group_by(unique_sim_number, frequency) %>%
  mutate(cumulative_incidence = cumsum(Incidence))

# Calculate mean cumulative incidence for each time and frequency
mean_cumulative_by_frequency <- analyticdf_cumulative %>%
  group_by(Time, frequency) %>%
  summarise(
    mean_cumulative_incidence = mean(cumulative_incidence, na.rm = TRUE),
    sd_cumulative_incidence = sd(cumulative_incidence, na.rm = TRUE),
    se_cumulative_incidence = sd(cumulative_incidence, na.rm = TRUE) / sqrt(n())
  )

# Plot the mean cumulative incidence over time by frequency
ggplot(mean_cumulative_by_frequency, aes(x = Time, y = mean_cumulative_incidence, color = as.factor(frequency))) +
  geom_line() +
  scale_x_continuous(breaks = seq(0,100,by = 10)) +
  labs(title = "Mean Cumulative Incidence Over Time by Frequency",
       x = "Time",
       y = "Mean Cumulative Incidence",
       color = "Frequency",
       fill = "Frequency") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_line(color = "black"),
    legend.position = "bottom"
  )

####Prop Outbreak50+ Over Time by Freq####
#Calculate the proportion of simulations with a cumulative 
#incidence of 50 or more at each time point for each frequency.

#Calculate Cumulative Incidence for Each Simulation and Time
analyticdf <- analyticdf %>%
  arrange(unique_sim_number, Time, frequency) %>%
  group_by(unique_sim_number, frequency) %>%
  mutate(cumulative_incidence = cumsum(Incidence))

#Calculate the Proportion of Simulations with Cumulative Incidence of 50 or More 
#at Each Time and Frequency
proportion_df <- analyticdf %>%
  group_by(Time, frequency) %>%
  summarise(
    num_sims = n_distinct(unique_sim_number),
    num_sims_50_or_more = sum(cumulative_incidence >= 50),
    proportion_50_or_more = num_sims_50_or_more / num_sims
  )
# Plot the proportion of simulations with cumulative incidence of 50 or more over time
ggplot(proportion_df, aes(x = Time, y = proportion_50_or_more, color = as.factor(frequency))) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(limits = c(0, 0.5)) +  # Since it's a proportion, y-axis ranges from 0 to 1
  labs(
    title = "Proportion of Simulations with Cumulative Incidence of 50 or More",
    x = "Time",
    y = "Proportion",
    color = "Frequency"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_line(color = "black"),
    legend.position = "bottom"
  )



####Proportion Zero Incidence by Frequency####
proportion_zero_df <- analyticdf %>%
  group_by(Time, frequency) %>%
  summarise(
    num_sims = n_distinct(unique_sim_number),
    num_sims_zero_incidence = sum(Incidence == 0),
    proportion_zero_incidence = num_sims_zero_incidence / num_sims
  )
# Plot the proportion of simulations with zero incidence over time
ggplot(proportion_zero_df, aes(x = Time, y = proportion_zero_incidence, color = as.factor(frequency))) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(limits = c(0, 1)) +  # Since it's a proportion, y-axis ranges from 0 to 1
  labs(
    title = "Proportion of Simulations with Zero Incidence",
    x = "Time",
    y = "Proportion",
    color = "Frequency"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_line(color = "black"),
    legend.position = "bottom"
  )



####Proportion Zero Incidence over time####
# Group the data by Time and calculate the proportion of simulations with zero Incidence
proportion_zero_incidence <- analyticdf %>%
  group_by(Time) %>%
  summarise(proportion_zero = mean(Incidence == 0))
# Plot the proportion of zero incidence over time
ggplot(proportion_zero_incidence, aes(x = Time, y = proportion_zero)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  labs(title = "Proportion of Simulations with Zero Incidence",
       x = "Time",
       y = "Proportion of Simulations with Zero Incidence") +
  theme_minimal()




#Instructions: https://bookdown.org/content/b298e479-b1ab-49fa-b83d-a57c2b034d49/correlation.html#using-the-levelplot-function-of-lattice 

####Create OutbreakProbs####
#for each frequency and test sensitivity combo, calculate the probability of an outbreak50
#tidyverse approach==> mutate with groups 

OutbreakProbs<-minianalyticdf %>% 
  group_by(frequency, sensitivity, isolationdelay) %>% summarize(proboutbreak50=mean(Outbreak50), proboutbreak20=mean(Outbreak20), proboutbreak10=mean(Outbreak10))

####Scenario Outbreak Probabilities####
NoRoutineTesting<-mean(OutbreakProbs$proboutbreak50[OutbreakProbs$frequency==0])
LowIntensityPessimistic<- OutbreakProbs$proboutbreak50[OutbreakProbs$frequency==1 & OutbreakProbs$sensitivity==50 & OutbreakProbs$isolationdelay==2]
LowIntensityOptimistic<- OutbreakProbs$proboutbreak50[OutbreakProbs$frequency==1 & OutbreakProbs$sensitivity==80 & OutbreakProbs$isolationdelay==1]
ModerateIntensityPessimistic<- OutbreakProbs$proboutbreak50[OutbreakProbs$frequency==2 & OutbreakProbs$sensitivity==50 & OutbreakProbs$isolationdelay==2]
ModerateIntensityOptimistic<- OutbreakProbs$proboutbreak50[OutbreakProbs$frequency==2 & OutbreakProbs$sensitivity==80 & OutbreakProbs$isolationdelay==1]
HighIntensityPessimistic<- OutbreakProbs$proboutbreak50[OutbreakProbs$frequency==5 & OutbreakProbs$sensitivity==50 & OutbreakProbs$isolationdelay==2]
HighIntensityOptimistic<- OutbreakProbs$proboutbreak50[OutbreakProbs$frequency==5 & OutbreakProbs$sensitivity==80 & OutbreakProbs$isolationdelay==1]
ScenarioOutbreakProbs<- as.data.frame(cbind(NoRoutineTesting, LowIntensityPessimistic, LowIntensityOptimistic, ModerateIntensityPessimistic, ModerateIntensityOptimistic, HighIntensityPessimistic, HighIntensityOptimistic))
rm(NoRoutineTesting, LowIntensityPessimistic, LowIntensityOptimistic, ModerateIntensityPessimistic, ModerateIntensityOptimistic, HighIntensityPessimistic, HighIntensityOptimistic)
ScenarioOutbreakProbs<- round(ScenarioOutbreakProbs,2)
write.csv(ScenarioOutbreakProbs, file = "ScenarioOutbreakProbs.csv")

#### Basic Heatmap 50#### 
basicheatmap50<-ggplot(OutbreakProbs, aes(frequency, sensitivity, fill= proboutbreak50)) + 
  geom_tile()+ scale_fill_gradient(low="blue", high="red")+ theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"))
basicheatmap50

####Fancy Heatmaps####
# alternative color scheme:https://r-graph-gallery.com/27-levelplot-with-lattice  
#levelplot(proboutbreak50 ~ frequency * sensitivity, OutbreakProbs, 
#          panel = panel.levelplot.points, cex = 0.01, col.regions = topo.colors(1000)
#    ) + 
#    layer_(panel.2dsmoother(..., n = 1000))

# showing data points on the same color scale 
x.scale <- list(at=seq(0,7,1))
y.scale <- list(at=seq(50,100,10))

col <- colorRampPalette(brewer.pal(8, "YlOrRd"))(25)

####Freq*sens Heatmap 50+ outbreak threshold####
# save the image you're about to create as a png image in specific directory with 600*350 resolution
png(file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/heatmapfreqsens50.png",
    width=600, height=350)
#create the image
freqsens50<-levelplot(proboutbreak50 ~ frequency * sensitivity, OutbreakProbs, 
                      panel = panel.levelplot.points, cex = 0.01, col.regions = col, xlab = "Test Frequency Per Week", ylab = "Test Sensitivity %", colorkey = TRUE, main= "Probability of an Outbreak of 50+", scales=list(x=x.scale, y=y.scale)
) + layer_(panel.2dsmoother(..., n = 1000))
freqsens50
# call dev.off() function to save the file
dev.off()

####Freq*sens Heatmap 20+ outbreak threshold####
# save the image you're about to create as a png image in specific directory with 600*350 resolution
png(file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/freqsens20.png",
    width=600, height=350)
#create the image
freqsens20<-levelplot(proboutbreak20 ~ frequency * sensitivity, OutbreakProbs, 
                      panel = panel.levelplot.points, cex = 0.01, col.regions = col, xlab = "Test Frequency Per Week", ylab = "Test Sensitivity %", colorkey = TRUE, main= "Probability of an Outbreak of 20+", scales=list(x=x.scale, y=y.scale)
) + layer_(panel.2dsmoother(..., n = 1000))
freqsens20
# call dev.off() function to save the file
dev.off()

####Freq*sens Heatmap 10+ outbreak threshold####
# save the image you're about to create as a png image in specific directory with 600*350 resolution
png(file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/freqsens10.png",
    width=600, height=350)
#create the image
freqsens10<-levelplot(proboutbreak10 ~ frequency * sensitivity, OutbreakProbs, 
                      panel = panel.levelplot.points, cex = 0.01, col.regions = col, xlab = "Test Frequency Per Week", ylab = "Test Sensitivity %", colorkey = TRUE, main= "Probability of an Outbreak of 10+", scales=list(x=x.scale, y=y.scale)
) + layer_(panel.2dsmoother(..., n = 1000))
freqsens10
# call dev.off() function to save the file
dev.off()

####Freq*delay Heatmap 50+ outbreak threshold####

# showing data points on the same color scale 
x.scale <- list(at=seq(0,7,1))
y.scale <- list(at=seq(0,3,1))

col <- colorRampPalette(brewer.pal(8, "YlOrRd"))(25)

# save the image you're about to create as a png image in specific directory with 600*350 resolution
png(file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/heatmapfreqdelay50.png",
    width=600, height=350)
#create the image
freqisodelay50<-levelplot(proboutbreak50 ~ frequency * isolationdelay, OutbreakProbs,
                          panel = panel.levelplot.points, cex = 0.01, col.regions = col, xlab = "Test Frequency Per Week", ylab = "Avg Length of Isolation Delay in Days", colorkey = TRUE, main= "Probability of an Outbreak of 50+", scales=list(x=x.scale, y=y.scale)
) + layer_(panel.2dsmoother(..., n = 1000))
freqisodelay50
# call dev.off() function to save the file
dev.off()


####Freq*delay Heatmap 20+ outbreak threshold####
# save the image you're about to create as a png image in specific directory with 600*350 resolution
png(file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/freqdelay20.png",
    width=600, height=350)
#create the image
freqisodelay20<-levelplot(proboutbreak20 ~ frequency * isolationdelay, OutbreakProbs,
                          panel = panel.levelplot.points, cex = 0.01, col.regions = col, xlab = "Test Frequency Per Week", ylab = "Avg Length of Isolation Delay in Days", colorkey = TRUE, main= "Probability of an Outbreak of 20+", scales=list(x=x.scale, y=y.scale)
) + layer_(panel.2dsmoother(..., n = 1000))
freqisodelay20
# call dev.off() function to save the file
dev.off()

####Freq*delay Heatmap 10+ outbreak threshold####
# pre-save the image you're about to create as a png image in specific directory with 600*350 resolution
png(file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/freqdelay10.png",
    width=600, height=350)
#create the image
freqisodelay10<-levelplot(proboutbreak10 ~ frequency * isolationdelay, OutbreakProbs,
                          panel = panel.levelplot.points, cex = 0.01, col.regions = col, xlab = "Test Frequency Per Week", ylab = "Avg Length of Isolation Delay in Days", colorkey = TRUE, main= "Probability of an Outbreak of 10+", scales=list(x=x.scale, y=y.scale)
) + layer_(panel.2dsmoother(..., n = 1000))
freqisodelay10
# call dev.off() function to save the file
dev.off()

####Create Long DF for Plots####
#reshape the data to long format
LongOutbreakProbs <- melt(OutbreakProbs, id = c("frequency","sensitivity","isolationdelay"), value.name = "OutbreakProbability")

#modify the dataframe so the variables are clear
# Rename the column and the values in the factor
levels(LongOutbreakProbs$variable)[levels(LongOutbreakProbs$variable)=="proboutbreak50"] <- "Probability of a 50+ outbreak"
levels(LongOutbreakProbs$variable)[levels(LongOutbreakProbs$variable)=="proboutbreak20"] <- "Probability of a 20+ outbreak"
levels(LongOutbreakProbs$variable)[levels(LongOutbreakProbs$variable)=="proboutbreak10"] <- "Probability of a 10+ outbreak"
names(LongOutbreakProbs)[names(LongOutbreakProbs)=="variable"]  <- "Threshold"

####FreqUnifiedMarginalPlots####
# pre-save the image you're about to create as a png image in specific directory with 600*350 resolution
png(file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/frequnifmarginalplots.png",
    width=600, height=350)

#create the plot
freq_unif_plots <- ggplot(LongOutbreakProbs,            
                          aes(x = frequency,
                              y = OutbreakProbability,
                              color = Threshold)) +
  labs(title = "Probability of an Outbreak by Test Frequency",
        x = "Test Frequency Per Week", y= "Probability of an Outbreak",caption = "\n*Minor horizontal jitter added for visibility, test frequency was always a multiple of 1.")+
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6,7,8)) +
  geom_smooth(se=FALSE)+   geom_jitter(width = 0.1, height=0)+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"))

freq_unif_plots
# call dev.off() function to save the file
dev.off()


####SensUnifiedMarginalPlots####
# pre-save the image you're about to create as a png image in specific directory with 600*350 resolution
png(file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/sensunifmarginalplots.png",
    width=600, height=350)
#create the plot
sens_unif_plots <- ggplot(LongOutbreakProbs,            
                          aes(x = sensitivity,
                              y = OutbreakProbability,
                              color = Threshold)) +
  labs(title = "Probability of an Outbreak by Test Sensitivity", x= "Test Sensitivity (%)", y= "Probability of an Outbreak", caption = "\n*Minor horizontal jitter added for visibility, test sensitivity was always a multiple of 10%.")+
  scale_x_continuous(breaks = c(50, 60, 70, 80, 90, 100)) + 
  scale_y_continuous(breaks= c(0,0.1,0.2,0.3,0.4, 0.5, 0.6)) +geom_smooth(se=FALSE) + 
  geom_jitter(width = 1, height=0)+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"))

sens_unif_plots
# call dev.off() function to save the file
dev.off()


####IsoDelayUnifiedMarginalPlots####
# pre-save the image you're about to create as a png image in specific directory with 600*350 resolution
png(file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/delayunifmarginalplots.png",
    width=600, height=350)
#create the plot
isodelay_unif_plots <- ggplot(LongOutbreakProbs,            
                              aes(x = isolationdelay,
                                  y = OutbreakProbability,
                                  color = Threshold))+
 labs(title = "Probability of an Outbreak by Isolation Delay", x= "Avg Isolation Delay in Days", y= "Probability of an Outbreak", caption = "\n*Minor horizontal jitter added for visibility, average isolation delays varied by 1.")+
scale_x_continuous(breaks = c(0, 1, 2, 3)) + geom_smooth(se=FALSE) + geom_jitter(width = 0.1, height=0)+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"))

isodelay_unif_plots
# call dev.off() function to save the file
dev.off()


####3D Plot####
#Temporarily edit names of variables for graphing
#colnames(OutbreakProbs)<- c("frequency", "sensitivity", "isolationdelay", "Pr(50+Infections)", "Pr(20+Infections)", "Pr(10+Infections)")
threeDplotA20 <- plot_ly(data = OutbreakProbs, 
                         x = ~frequency, 
                         y = ~sensitivity, 
                         z = ~isolationdelay) %>%
  add_markers(
    color = ~proboutbreak20,
    marker = list(colorbar = list(title = ""))
  ) %>%
  layout(
    title = "Probability of Outbreak of 20+",
    scene = list(
      xaxis = list(title = "Tests Per Week"),
      yaxis = list(title = "Test Sensitivity (%)"),
      zaxis = list(title = "Avg Isolation Delay")
    )
  )

threeDplotA20

####**CREATE TABLES**####

####Modify vars for logit regression####
#center sensitivity at its mean
minianalyticdf$sensitivity<- minianalyticdf$sensitivity-mean(minianalyticdf$sensitivity) 


####Simple Logit50####
simplelogit50<-glm(Outbreak50~frequency+sensitivity+isolationdelay, data=minianalyticdf, family=binomial(link = "logit"))
summary(simplelogit50)
OddsRatios<-exp(coef(simplelogit50))
names(OddsRatios)<- c("Intercept","Frequency","Sensitivity","Isolation Delay")
OddsRatios<- round(OddsRatios,3)
Sumsimplelogit50<- summary(simplelogit50)
logOR<-Sumsimplelogit50$coefficients[,1]
SEs<- Sumsimplelogit50$coefficients[,2]
Z_score<- Sumsimplelogit50$coefficients[,3]
P_values<-Sumsimplelogit50$coefficients[,4]
simplelogit50table<- as.data.frame(cbind(OddsRatios,logOR,SEs,Z_score, P_values))
simplelogit50table$LowerBoundCI<-exp(logOR-1.96*SEs)
simplelogit50table$UpperBoundCI<-exp(logOR+1.96*SEs)
simplelogit50table<- simplelogit50table[,c(1,6,7,4,5)]
simplelogit50table[,1:4] <- round(simplelogit50table[,1:4], 3)
simplelogit50table[,5]<- "<0.001"
names(simplelogit50table)<- c("Odds Ratios", "Lower Bound 95% CI", "Upper Bound 95% CI", "Z-Score", "P-Value")
write.csv(simplelogit50table, file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/simplelogit50table.csv", row.names=TRUE)

####Logit With Interaction 50####
logitwint50<-glm(Outbreak50~frequency + sensitivity + isolationdelay + frequency*sensitivity + frequency*isolationdelay + sensitivity*isolationdelay, data=minianalyticdf, family=binomial(link = "logit"))
summary(logitwint50)
OddsRatios2<-exp(coef(logitwint50))
names(OddsRatios2)<- c("Intercept","Frequency","Sensitivity","Isolation Delay","Frequency:Sensitivity","Frequency:Isolation Delay", "Sensitivity:Isolation Delay")
OddsRatios2<- round(OddsRatios2,3)
Sumlogitwint50<- summary(logitwint50)
logOR2<-Sumlogitwint50$coefficients[,1]
SEs2<- Sumlogitwint50$coefficients[,2]
Z_score2<- Sumlogitwint50$coefficients[,3]
P_values2<-Sumlogitwint50$coefficients[,4]
logitwint50table<- as.data.frame(cbind(OddsRatios2,logOR2,SEs2,Z_score2, P_values2))
logitwint50table$LowerBoundCI<-exp(logOR2-1.96*SEs2)
logitwint50table$UpperBoundCI<-exp(logOR2+1.96*SEs2)
logitwint50table<- logitwint50table[,c(1,6,7,4,5)]
logitwint50table[,1:4]<- round(logitwint50table,3)
logitwint50table[,5]<- "<0.001"
names(logitwint50table)<- c("Odds Ratios", "Lower Bound 95% CI", "Upper Bound 95% CI", "Z-Score", "P-Value")
write.csv(logitwint50table, file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/logitwint50table.csv", row.names=TRUE)

####Table ORFreqbyIsoDelay####

# Create a new dataframe to store results
Tests_Per_Week <- 0:7
ORforFrequencybyIsoDelay <- data.frame(Tests_Per_Week)

# Function to calculate values for the table
calc_values <- function(or_freq, lb_freq, ub_freq, or_fid, lb_fid, ub_fid, t) {
  or_0d <- or_freq^t
  lb_0d <- lb_freq^t
  ub_0d <- ub_freq^t
  
  or_1d <- (or_freq * or_fid)^t
  lb_1d <- (lb_freq * lb_fid)^t
  ub_1d <- (ub_freq * ub_fid)^t
  
  or_2d <- (or_freq * (or_fid^2))^t
  lb_2d <- (lb_freq * (lb_fid^2))^t
  ub_2d <- (ub_freq * (ub_fid^2))^t
  
  or_3d <- (or_freq * (or_fid^3))^t
  lb_3d <- (lb_freq * (lb_fid^3))^t
  ub_3d <- (ub_freq * (ub_fid^3))^t
  
  return(c(or_0d, lb_0d, ub_0d, or_1d, lb_1d, ub_1d, or_2d, lb_2d, ub_2d, or_3d, lb_3d, ub_3d))
}

# Extract the relevant values from logitwint50table
or_freq <- logitwint50table["Frequency", "Odds Ratios"]
lb_freq <- logitwint50table["Frequency", "Lower Bound 95% CI"]
ub_freq <- logitwint50table["Frequency", "Upper Bound 95% CI"]

or_fid <- logitwint50table["Frequency:Isolation Delay", "Odds Ratios"]
lb_fid <- logitwint50table["Frequency:Isolation Delay", "Lower Bound 95% CI"]
ub_fid <- logitwint50table["Frequency:Isolation Delay", "Upper Bound 95% CI"]

# Calculate values for each row
for (i in 1:nrow(ORforFrequencybyIsoDelay)) {
  t <- ORforFrequencybyIsoDelay$Tests_Per_Week[i]
  ORforFrequencybyIsoDelay[i, 2:13] <- calc_values(or_freq, lb_freq, ub_freq, or_fid, lb_fid, ub_fid, t)
}

# Set column names
colnames(ORforFrequencybyIsoDelay) <- c(
  "Tests Per Week",
  "OR with Avg Isolation Delay of 0 Days", "Lower Bound 95% CI OR with Avg Isolation Delay of 0 Days", "Upper Bound 95% CI OR with Avg Isolation Delay of 0 Days",
  "OR with Avg Isolation Delay of 1 Day", "Lower Bound 95% CI OR with Avg Isolation Delay of 1 Day", "Upper Bound 95% CI OR with Avg Isolation Delay of 1 Day",
  "OR with Avg Isolation Delay of 2 Days", "Lower Bound 95% CI OR with Avg Isolation Delay of 2 Days", "Upper Bound 95% CI OR with Avg Isolation Delay of 2 Days",
  "OR with Avg Isolation Delay of 3 Days", "Lower Bound 95% CI OR with Avg Isolation Delay of 3 Days", "Upper Bound 95% CI OR with Avg Isolation Delay of 3 Days"
)

# Display the new table
print(ORforFrequencybyIsoDelay)

#round it to 3 digits
ORforFrequencybyIsoDelay<- round(ORforFrequencybyIsoDelay,3)
#transpose it
#ORforFrequencybyIsoDelay<- t(ORforFrequencybyIsoDelay)
#save it
#write.csv(ORforFrequencybyIsoDelay, file = "ORforFrequencybyIsoDelay.csv")

### next step: combining columns (unsucessful-hence commented out)###

# # Function to combine OR and CI values into the desired format
# combine_or_ci <- function(or, lb, ub) {
#   paste0(or, " (", lb, " to ", ub, ")")
# }
# 
# # Create a new dataframe to store the combined results
# TableORforFrequencybyIsoDelay <- data.frame(
#   Tests_Per_Week = ORforFrequencybyIsoDelay$Tests_Per_Week,
#   `OR with Avg Isolation Delay of 0 Days` = character(nrow(ORforFrequencybyIsoDelay)),
#   `OR with Avg Isolation Delay of 1 Day` = character(nrow(ORforFrequencybyIsoDelay)),
#   `OR with Avg Isolation Delay of 2 Days` = character(nrow(ORforFrequencybyIsoDelay)),
#   `OR with Avg Isolation Delay of 3 Days` = character(nrow(ORforFrequencybyIsoDelay)),
#   stringsAsFactors = FALSE
# )
# 
# # Combine the values for each column set
# for (i in 1:nrow(ORforFrequencybyIsoDelay)) {
#   TableORforFrequencybyIsoDelay$`OR with Avg Isolation Delay of 0 Days`[i] <- combine_or_ci(
#     ORforFrequencybyIsoDelay[i, 2], ORforFrequencybyIsoDelay[i, 3], ORforFrequencybyIsoDelay[i, 4])
#   
#   TableORforFrequencybyIsoDelay$`OR with Avg Isolation Delay of 1 Day`[i] <- combine_or_ci(
#     ORforFrequencybyIsoDelay[i, 5], ORforFrequencybyIsoDelay[i, 6], ORforFrequencybyIsoDelay[i, 7])
#   
#   TableORforFrequencybyIsoDelay$`OR with Avg Isolation Delay of 2 Days`[i] <- combine_or_ci(
#     ORforFrequencybyIsoDelay[i, 8], ORforFrequencybyIsoDelay[i, 9], ORforFrequencybyIsoDelay[i, 10])
#   
#   TableORforFrequencybyIsoDelay$`OR with Avg Isolation Delay of 3 Days`[i] <- combine_or_ci(
#     ORforFrequencybyIsoDelay[i, 11], ORforFrequencybyIsoDelay[i, 12], ORforFrequencybyIsoDelay[i, 13])
# }
# 
# # Display the new table to verify
# print(TableORforFrequencybyIsoDelay)




####Graph&Table Freq:Delay50####
freqbydelay50int<- as.data.frame(cbind(frequency=(c(seq(0,7,1))), 
                                  ORdelay0=c(1,logitwint50table[2,1], logitwint50table[2,1]^2, logitwint50table[2,1]^3, logitwint50table[2,1]^4, logitwint50table[2,1]^5, logitwint50table[2,1]^6, logitwint50table[2,1]^7), 
                                  ORdelay1= c(1, logitwint50table[2,1]*logitwint50table[6,1],logitwint50table[2,1]^2*logitwint50table[6,1], logitwint50table[2,1]^3*logitwint50table[6,1], logitwint50table[2,1]^4*logitwint50table[6,1], logitwint50table[2,1]^5*logitwint50table[6,1], logitwint50table[2,1]^6*logitwint50table[6,1], logitwint50table[2,1]^7*logitwint50table[6,1]),
                                  ORdelay2= c(1, logitwint50table[2,1]*logitwint50table[6,1]^2,logitwint50table[2,1]^2*logitwint50table[6,1]^2, logitwint50table[2,1]^3*logitwint50table[6,1]^2, logitwint50table[2,1]^4*logitwint50table[6,1]^2, logitwint50table[2,1]^5*logitwint50table[6,1]^2, logitwint50table[2,1]^6*logitwint50table[6,1]^2, logitwint50table[2,1]^7*logitwint50table[6,1]^2),
                                  ORdelay3= c(1, logitwint50table[2,1]*logitwint50table[6,1]^3,logitwint50table[2,1]^2*logitwint50table[6,1]^3, logitwint50table[2,1]^3*logitwint50table[6,1]^3, logitwint50table[2,1]^4*logitwint50table[6,1]^3, logitwint50table[2,1]^5*logitwint50table[6,1]^3, logitwint50table[2,1]^6*logitwint50table[6,1]^3, logitwint50table[2,1]^7*logitwint50table[6,1]^3)
))
long_freqbydelay50int <- melt(freqbydelay50int, id.vars = "frequency", variable.name = "Isolation_Delay", value.name = "OR")

# Map the labels for Isolation_Delay
long_freqbydelay50int$Isolation_Delay <- factor(long_freqbydelay50int$Isolation_Delay,
                                           levels = c("ORdelay0","ORdelay1","ORdelay2", "ORdelay3"),
                                           labels = c("Isolation delay of 0 days","Isolation delay of 1 day", "Isolation delay of 2 days", "Isolation delay of 3 days"))

#clean up freqbydelay50int table
names(freqbydelay50int)<- c("Tests Per Week", "OR with Avg Isolation Delay of 0 Days", "OR with Avg Isolation Delay of 1 Day", "OR with Avg Isolation Delay of 2 Days", "OR with Avg Isolation Delay of 3 Days")
freqbydelay50int<- round(freqbydelay50int,2)
#save freqbydelay50int table
write.csv(freqbydelay50int, file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/freqbydelay50int.csv", row.names = FALSE)

# Create plotfreqbydelay50int 

# save the image you're about to create as a png image in specific directory with 600*350 resolution
png(file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/plotfreqbydelay50int.png",
    width=600, height=350)

plotfreqbydelay50int<-ggplot(long_freqbydelay50int, aes(x = frequency, y = OR, color = Isolation_Delay, linetype = Isolation_Delay)) +
  geom_point() +
  geom_line(size = 1.5) +
  labs(title = "Isolation Delay Moderates Test Frequency",
       x = "Test Frequency Per Week",
       y = "Odds Ratio (OR) for Outbreak of 50+",
       color = "Isolation Delay",
       linetype = "Isolation Delay") +
  scale_x_continuous(breaks = seq(0, 7, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_linetype_manual(values = rep("solid", length(unique(long_freqbydelay50int$Isolation_Delay)))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0, size = 22, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"))
plotfreqbydelay50int

# call dev.off() function to save the file
dev.off()

####Table Sens:Delay####

sensbydelay50int<- as.data.frame(cbind(ORdelay0=c(logitwint50table[3,1]), 
                                       ORdelay1= c(logitwint50table[3,1]*logitwint50table[7,1]),
                                       ORdelay2= c(logitwint50table[3,1]*logitwint50table[7,1]^2),
                                       ORdelay3= c(logitwint50table[3,1]*logitwint50table[7,1]^3)
))
names(sensbydelay50int)<- c("OR with Avg Isolation Delay of 0 Days", "OR with Avg Isolation Delay of 1 Day", "OR with Avg Isolation Delay of 2 Days", "OR with Avg Isolation Delay of 3 Days")
sensbydelay50int<-round(sensbydelay50int,3)
write.csv(sensbydelay50int,file= "/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/sensbydelay50int.csv", row.names = FALSE )




####GraphSens:Freq#### check again
sensbyfreq50int<- as.data.frame(cbind(sensitivity=(c(seq(50,100,10), 
                                      ORfreq0=c(1, logitwint50table[3,1]^10, logitwint50table[3,1]^10^2, logitwint50table[3,1]^10^3, logitwint50table[3,1]^10^4, logitwint50table[3,1]^10^5, logitwint50table[3,1]^10^6, logitwint50table[3,1]^10^7), 
                                      ORfreq1= c(1, logitwint50table[3,1]^10*logitwint50table[6,1],logitwint50table[3,1]^10^2*logitwint50table[6,1], logitwint50table[3,1]^10^3*logitwint50table[6,1], logitwint50table[3,1]^10^4*logitwint50table[6,1], logitwint50table[3,1]^10^5*logitwint50table[6,1], logitwint50table[3,1]^10^6*logitwint50table[6,1], logitwint50table[3,1]^10^7*logitwint50table[6,1]),
                                      ORfreq2= c(1, logitwint50table[3,1]^10*logitwint50table[6,1]^2,logitwint50table[3,1]^10^2*logitwint50table[6,1]^2, logitwint50table[3,1]^10^3*logitwint50table[6,1]^2, logitwint50table[3,1]^10^4*logitwint50table[6,1]^2, logitwint50table[3,1]^10^5*logitwint50table[6,1]^2, logitwint50table[3,1]^10^6*logitwint50table[6,1]^2, logitwint50table[3,1]^10^7*logitwint50table[6,1]^2),
                                      ORfreq3= c(1, logitwint50table[3,1]^10*logitwint50table[6,1]^3,logitwint50table[3,1]^10^2*logitwint50table[6,1]^3, logitwint50table[3,1]^10^3*logitwint50table[6,1]^3, logitwint50table[3,1]^10^4*logitwint50table[6,1]^3, logitwint50table[3,1]^10^5*logitwint50table[6,1]^3, logitwint50table[3,1]^10^6*logitwint50table[6,1]^3, logitwint50table[3,1]^10^7*logitwint50table[6,1]^3),
                                      ORfreq4= c(1, logitwint50table[3,1]^10*logitwint50table[6,1]^4,logitwint50table[3,1]^10^2*logitwint50table[6,1]^4, logitwint50table[3,1]^10^3*logitwint50table[6,1]^4, logitwint50table[3,1]^10^4*logitwint50table[6,1]^4, logitwint50table[3,1]^10^5*logitwint50table[6,1]^4, logitwint50table[3,1]^10^6*logitwint50table[6,1]^4, logitwint50table[3,1]^10^7*logitwint50table[6,1]^4),
                                      ORfreq5= c(1, logitwint50table[3,1]^10*logitwint50table[6,1]^5,logitwint50table[3,1]^10^2*logitwint50table[6,1]^5, logitwint50table[3,1]^10^3*logitwint50table[6,1]^5, logitwint50table[3,1]^10^4*logitwint50table[6,1]^5, logitwint50table[3,1]^10^5*logitwint50table[6,1]^5, logitwint50table[3,1]^10^6*logitwint50table[6,1]^5, logitwint50table[3,1]^10^7*logitwint50table[6,1]^5),
                                      ORfreq6= c(1, logitwint50table[3,1]^10*logitwint50table[6,1]^6,logitwint50table[3,1]^10^2*logitwint50table[6,1]^6, logitwint50table[3,1]^10^3*logitwint50table[6,1]^6, logitwint50table[3,1]^10^4*logitwint50table[6,1]^6, logitwint50table[3,1]^10^5*logitwint50table[6,1]^6, logitwint50table[3,1]^10^6*logitwint50table[6,1]^6, logitwint50table[3,1]^10^7*logitwint50table[6,1]^6),
                                      ORfreq7= c(1, logitwint50table[3,1]^10*logitwint50table[6,1]^7,logitwint50table[3,1]^10^2*logitwint50table[6,1]^7, logitwint50table[3,1]^10^3*logitwint50table[6,1]^7, logitwint50table[3,1]^10^4*logitwint50table[6,1]^7, logitwint50table[3,1]^10^5*logitwint50table[6,1]^7, logitwint50table[3,1]^10^6*logitwint50table[6,1]^7, logitwint50table[3,1]^10^7*logitwint50table[6,1]^7),
                                                                            ))))

long_sensbyfreq50int <- melt(sensbyfreq50int, id.vars = "sensitivity", variable.name = "Test Frequency Per Week", value.name = "OR")

# Map the labels for Isolation_Delay
long_freqbydelay50int$Isolation_Delay <- factor(long_freqbydelay50int$Isolation_Delay,
                                                levels = c("ORdelay0","ORdelay1","ORdelay2", "ORdelay3"),
                                                labels = c("Isolation delay of 0 days","Isolation delay of 1 day", "Isolation delay of 2 days", "Isolation delay of 3 days"))


####Simple Logit20####
simplelogit20<-glm(Outbreak20~frequency+sensitivity+isolationdelay, data=minianalyticdf, family=binomial(link = "logit"))
summary(simplelogit20)
OddsRatios<-exp(coef(simplelogit20))
names(OddsRatios)<- c("Intercept","Frequency","Sensitivity","Isolation Delay")
OddsRatios<- round(OddsRatios,3)
Sumsimplelogit20<- summary(simplelogit20)
logOR<-Sumsimplelogit20$coefficients[,1]
SEs<- Sumsimplelogit20$coefficients[,2]
Z_score<- Sumsimplelogit20$coefficients[,3]
P_values<-Sumsimplelogit20$coefficients[,4]
simplelogit20table<- as.data.frame(cbind(OddsRatios,logOR,SEs,Z_score, P_values))
simplelogit20table$LowerBoundCI<-exp(logOR-1.96*SEs)
simplelogit20table$UpperBoundCI<-exp(logOR+1.96*SEs)
simplelogit20table<- simplelogit20table[,c(1,6,7,4,5)]
simplelogit20table<- round(simplelogit20table,3)
names(simplelogit20table)<- c("Odds Ratios", "Lower Bound 95% CI", "Upper Bound 95% CI", "Z-Score", "P-Value")
write.csv(simplelogit20table, file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/simplelogit20table.csv", row.names=TRUE)

####Logit With Interaction 20####
logitwint20<-glm(Outbreak20~frequency + sensitivity + isolationdelay + frequency*sensitivity + frequency*isolationdelay + sensitivity*isolationdelay, data=minianalyticdf, family=binomial(link = "logit"))
summary(logitwint20)
OddsRatios2<-exp(coef(logitwint20))
names(OddsRatios2)<- c("Intercept","Frequency","Sensitivity","Isolation Delay","Frequency:Sensitivity","Frequency:Isolation Delay", "Sensitivity:Isolation Delay")
OddsRatios2<- round(OddsRatios2,3)
Sumlogitwint20<- summary(logitwint20)
logOR2<-Sumlogitwint20$coefficients[,1]
SEs2<- Sumlogitwint20$coefficients[,2]
Z_score2<- Sumlogitwint20$coefficients[,3]
P_values2<-Sumlogitwint20$coefficients[,4]
logitwint20table<- as.data.frame(cbind(OddsRatios2,logOR2,SEs2,Z_score2, P_values2))
logitwint20table$LowerBoundCI<-exp(logOR2-1.96*SEs2)
logitwint20table$UpperBoundCI<-exp(logOR2+1.96*SEs2)
logitwint20table<- logitwint20table[,c(1,6,7,4,5)]
logitwint20table<- round(logitwint20table,3)
names(logitwint20table)<- c("Odds Ratios", "Lower Bound 95% CI", "Upper Bound 95% CI", "Z-Score", "P-Value")
write.csv(logitwint20table, file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/logitwint20table.csv", row.names=TRUE)

####Simple Logit10####
simplelogit10<-glm(Outbreak10~frequency+sensitivity+isolationdelay, data=minianalyticdf, family=binomial(link = "logit"))
summary(simplelogit10)
OddsRatios<-exp(coef(simplelogit10))
names(OddsRatios)<- c("Intercept","Frequency","Sensitivity","Isolation Delay")
OddsRatios<- round(OddsRatios,3)
Sumsimplelogit10<- summary(simplelogit10)
logOR<-Sumsimplelogit10$coefficients[,1]
SEs<- Sumsimplelogit10$coefficients[,2]
Z_score<- Sumsimplelogit10$coefficients[,3]
P_values<-Sumsimplelogit10$coefficients[,4]
simplelogit10table<- as.data.frame(cbind(OddsRatios,logOR,SEs,Z_score, P_values))
simplelogit10table$LowerBoundCI<-exp(logOR-1.96*SEs)
simplelogit10table$UpperBoundCI<-exp(logOR+1.96*SEs)
simplelogit10table<- simplelogit10table[,c(1,6,7,4,5)]
simplelogit10table<- round(simplelogit10table,3)
names(simplelogit10table)<- c("Odds Ratios", "Lower Bound 95% CI", "Upper Bound 95% CI", "Z-Score", "P-Value")
write.csv(simplelogit10table, file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/simplelogit10table.csv", row.names=TRUE)

####Logit With Interaction 10####
logitwint10<-glm(Outbreak10~frequency + sensitivity + isolationdelay + frequency*sensitivity + frequency*isolationdelay + sensitivity*isolationdelay, data=minianalyticdf, family=binomial(link = "logit"))
summary(logitwint10)
OddsRatios2<-exp(coef(logitwint10))
names(OddsRatios2)<- c("Intercept","Frequency","Sensitivity","Isolation Delay","Frequency:Sensitivity","Frequency:Isolation Delay", "Sensitivity:Isolation Delay")
OddsRatios2<- round(OddsRatios2,3)
Sumlogitwint10<- summary(logitwint10)
logOR2<-Sumlogitwint10$coefficients[,1]
SEs2<- Sumlogitwint10$coefficients[,2]
Z_score2<- Sumlogitwint10$coefficients[,3]
P_values2<-Sumlogitwint10$coefficients[,4]
logitwint10table<- as.data.frame(cbind(OddsRatios2,logOR2,SEs2,Z_score2, P_values2))
logitwint10table$LowerBoundCI<-exp(logOR2-1.96*SEs2)
logitwint10table$UpperBoundCI<-exp(logOR2+1.96*SEs2)
logitwint10table<- logitwint10table[,c(1,6,7,4,5)]
logitwint10table<- round(logitwint10table,3)
names(logitwint10table)<- c("Odds Ratios", "Lower Bound 95% CI", "Upper Bound 95% CI", "Z-Score", "P-Value")
write.csv(logitwint10table, file="/Users/byroncohen/Dropbox/HSPH PhD/Dissertation/ID Surveillance/Biosafety/Stochastic Network Model Code/FASRC_Code/PaperOne_4_2_2024/logitwint10table.csv", row.names=TRUE)



####Create Poisson Model for Outbreak Size####

#simpleoutbreaksize<-lm(everinfected~frequency+sensitivity+isolationdelay, data= minianalyticdf)
#summary(simpleoutbreaksize)

