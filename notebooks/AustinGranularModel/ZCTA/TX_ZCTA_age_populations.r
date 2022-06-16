# Austin, TX ZCTAs 2018 and 2019

rm(list=ls())
library(tidyverse)
library(tidycensus)
library(sf)

# load environment with US Census Bureau API key
# key request documentation available at https://www.census.gov/data/developers.html
load('/Users/kpierce/uscb_query_env.RData')

API_KEY <- Sys.getenv("CENSUS_API_KEY")
census_api_key(API_KEY)

load_all_variables <- function(year){
  
  # get all the variable data
  basic <- load_variables(year, 'acs5', cache=TRUE)
  basic$table <- 'basic'
  table_name <- sapply(basic$name, strsplit, split="_")
  
  data_profiles <- load_variables(year, 'acs5/profile', cache=TRUE)
  data_profiles$table <- 'data_profiles'
  
  survey <- load_variables(year, 'acs5/subject', cache=TRUE)
  survey$table <- 'survey'
  
  all <- rbind(basic, data_profiles)
  all <- rbind(all, survey)
  
  # remove extraneous colons in the label field
  all$label <- sapply(all$label, gsub, pattern=':', replacement='')
  
  # separate aggregate columns
  all <- separate(all, name, into=c('table_name'), sep='_', remove=FALSE)
  all$year <- year
  
  return(all)
}

### 2019 ###

v19 <- load_all_variables(2019)
v19_male_ages <- v19 %>% filter(grepl('Estimate!!Total!!Male!!', label) & concept=='SEX BY AGE')
v19_female_ages <- v19 %>% filter(grepl('Estimate!!Total!!Female!!', label) & concept=='SEX BY AGE')
age_vars <- c(v19_male_ages$name, v19_female_ages$name)

age_var_names <- bind_rows(v19_male_ages, v19_female_ages)

tx_age_by_zcta <- get_acs(geography='zcta', year=2019, state='TX', variables=age_vars, survey="acs5", show_call=TRUE, geometry=FALSE)
tx_age_by_zcta_labeled <- left_join(tx_age_by_zcta, age_var_names, by=c('variable'='name'))

tx_age_by_zcta_labeled <- separate(tx_age_by_zcta_labeled, label, into=c('x', 'y', 'gender', 'age'), sep='!!', remove=FALSE)
tx_age_by_zcta_labeled <- tx_age_by_zcta_labeled %>% select('GEOID', 'estimate', 'age')
tx_age_by_zcta_labeled <- tx_age_by_zcta_labeled %>% group_by(GEOID, age) %>% summarise(estimate = sum(estimate))

write.csv(tx_age_by_zcta_labeled, '~/epimodels/notebooks/AustinGranularModel/zcta/2019_TX_zcta_age_populations.csv')

### 2018 ###

v18 <- load_all_variables(2018)
v18_male_ages <- v18 %>% filter(grepl('Estimate!!Total!!Male!!', label) & concept=='SEX BY AGE')
v18_female_ages <- v18 %>% filter(grepl('Estimate!!Total!!Female!!', label) & concept=='SEX BY AGE')
age_vars_2018 <- c(v18_male_ages$name, v18_female_ages$name)

age_var_names_2018 <- bind_rows(v18_male_ages, v18_female_ages)

tx_age_by_zcta_2018 <- get_acs(geography='zcta', year=2018, state='TX', variables=age_vars, survey="acs5", show_call=TRUE, geometry=FALSE)
tx_age_by_zcta_2018_labeled <- left_join(tx_age_by_zcta_2018, age_var_names_2018, by=c('variable'='name'))

tx_age_by_zcta_2018_labeled <- separate(tx_age_by_zcta_2018_labeled, label, into=c('x', 'y', 'gender', 'age'), sep='!!', remove=FALSE)
tx_age_by_zcta_2018_labeled <- tx_age_by_zcta_2018_labeled %>% select('GEOID', 'estimate', 'age')
tx_age_by_zcta_2018_labeled <- tx_age_by_zcta_2018_labeled %>% group_by(GEOID, age) %>% summarise(estimate = sum(estimate))

write.csv(tx_age_by_zcta_2018_labeled, '~/epimodels/notebooks/AustinGranularModel/zcta/2018_TX_zcta_age_populations.csv')

### ZCTA GEOMETRIES ###

library(tigris)
zcta_2018 <- zctas(starts_with=c('75', '75', '77', '78', '78'), year=2018)
st_write(zcta_2018, '~/epimodels/notebooks/AustinGranularModel/zcta/2018_TX_zcta.shp')

zcta_2019 <- zctas(starts_with=c('75', '75', '77', '78', '78'), year=2019)
st_write(zcta_2019, '~/epimodels/notebooks/AustinGranularModel/zcta/2019_TX_zcta.shp')

