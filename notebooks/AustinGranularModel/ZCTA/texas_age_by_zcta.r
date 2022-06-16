rm(list=ls())
library(tidycensus)
library(tidyverse)

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
  #all <- separate(all, label, into=c('l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'l8'), sep='!!', remove=FALSE)
  all <- separate(all, name, into=c('table_name'), sep='_', remove=FALSE)
  all$year <- year

  return(all)
}

v19 <- load_all_variables(2019)

v19_male_ages <- v19 %>% filter(grepl('Estimate!!Total!!Male!!', label) & concept=='SEX BY AGE')
v19_female_ages <- v19 %>% filter(grepl('Estimate!!Total!!Female!!', label) & concept=='SEX BY AGE')
age_vars <- c(v19_male_ages$name, v19_female_ages$name)

age_var_names <- bind_rows(v19_male_ages, v19_female_ages)

API_KEY <- read.table('~/CooksProTX/us_census_api_key.txt', header=TRUE, colClasses='character')
census_api_key(API_KEY)

tx_age_by_zcta <- get_acs('zcta', year=2019, state='TX', variables=age_vars, survey="acs5", show_call=TRUE)
tx_age_by_zcta_labeled <- left_join(tx_age_by_zcta, age_var_names, by=c('variable'='name'))

tx_age_by_zcta_labeled <- separate(tx_age_by_zcta_labeled, label, into=c('x', 'y', 'gender', 'age'), sep='!!', remove=FALSE)
tx_age_by_zcta_labeled <- tx_age_by_zcta_labeled %>% select('GEOID', 'estimate', 'age')
tx_age_by_zcta_labeled <- tx_age_by_zcta_labeled %>% group_by(GEOID, age) %>% summarise(estimate = sum(estimate))

write.csv(tx_age_by_zcta_labeled, '~/COVID19/SchoolCatchmentDemo/2019_TX_ZCTA_age_populations.csv')

