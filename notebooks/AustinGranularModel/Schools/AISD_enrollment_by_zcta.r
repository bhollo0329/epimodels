# school district populations
rm(list=ls())
library(tidyverse)
library(tidycensus)
library(sf)

### 2019 ###

API_KEY <- Sys.getenv("CENSUS_API_KEY")
census_api_key(API_KEY)

# download the school district shapes, and extract the AISD shape to use as a filter (boundary is unchanged from 2018-19)
school_distr_shape <- get_acs(
  geography='school district (unified)',
  variables=c('B09001_001'),
  state='TX',
  year=2019,
  survey='acs5',
  show_call=TRUE,
  geometry=TRUE
)
aisd_shape = school_distr_shape %>% filter(NAME=='Austin Independent School District, Texas')

# load the ZCTA shapefiles for 2018 and 2019
tx_zcta_2018 <- st_read('~/epimodels/notebooks/AustinGranularModel/ZCTA/2018_TX_zcta.shp')
tx_zcta_2019 <- st_read('~/epimodels/notebooks/AustinGranularModel/ZCTA/2019_TX_zcta.shp')

# identify zip codes that intersect AISD
aisd_zcta_2018 <- st_intersection(aisd_shape, tx_zcta_2018)
aisd_zcta_2019 <- st_intersection(aisd_shape, tx_zcta_2019)

# collect the enrollment variables and download enrollment data for AISD ZCTAs
v18 <- load_variables(2018, dataset='acs5/subject')
school_enrollment_2018 <- v18 %>% filter(grepl('ENROLL', v18$concept)) %>% filter(grepl('Percent', label)) # %>% filter(!('graduate' %in% 'label'))

enrollment_2018 <- get_acs(
  geography='zcta',
  variables=school_enrollment_2018$name,
  state='TX',
  zcta=aisd_zcta_2018$GEOID10,
  year=2018,
  survey='acs5',
  show_call=TRUE,
  geometry=FALSE
)

# join variables with their names, pivot, and save
enrollment_2018_named <- left_join(enrollment_2018, school_enrollment_2018, by=c('variable'='name'))
enrollment_2018_named <- enrollment_2018_named %>% select(-c(moe, variable))

private_enrollment <- enrollment_2018_named %>% filter(grepl('private', label))
public_enrollment <- enrollment_2018_named %>% filter(grepl('public', label))

public_enrollment_wide <- public_enrollment %>% 
  pivot_wider(names_from='label', values_from=c('estimate'))

private_enrollment_wide <- private_enrollment %>% 
  pivot_wider(names_from='label', values_from=c('estimate'))

write.csv(
  public_enrollment_wide,
  '~/epimodels/notebooks/AustinGranularModel/Schools/data/2018_AISD_Public_School_Enrollment_Estimates_by_ZCTA.csv'
  )

write.csv(
  private_enrollment_wide,
  '~/epimodels/notebooks/AustinGranularModel/Schools/data/2018_AISD_Private_School_Enrollment_Estimates_by_ZCTA.csv'
)

# repeat for 2019
# collect the enrollment variables and download enrollment data for AISD ZCTAs
v19 <- load_variables(2019, dataset='acs5/subject')
school_enrollment_2019 <- v19 %>% filter(grepl('ENROLL', v19$concept)) %>% filter(grepl('Percent', label)) # %>% filter(!('graduate' %in% 'label'))

# geography hierarchies changed in 2019 so we have to retrieve data for all the ZCTAs...
enrollment_2019 <- get_acs(
  geography='zcta',
  variables=school_enrollment_2019$name,
  #state='TX',
  #zcta=aisd_zcta_2019$GEOID10,
  year=2019,
  survey='acs5',
  show_call=TRUE,
  geometry=FALSE
)
enrollment_2019 <- enrollment_2019 %>% separate(NAME, into=c('geograpy', 'GEOID10'), sep=' ', remove=FALSE)
enrollment_2019_filtered <- enrollment_2019 %>% filter(GEOID10 %in% aisd_zcta_2019$GEOID10)

# join variables with their names, pivot, and save
enrollment_2019_named <- left_join(enrollment_2019_filtered, school_enrollment_2019, by=c('variable'='name'))
enrollment_2019_named <- enrollment_2019_named %>% select(-c(moe, variable))

private_enrollment_2019 <- enrollment_2019_named %>% filter(grepl('private', label))
public_enrollment_2019 <- enrollment_2019_named %>% filter(grepl('public', label))

public_enrollment_2019_wide <- public_enrollment_2019 %>% 
  pivot_wider(names_from='label', values_from=c('estimate'))

private_enrollment_2019_wide <- private_enrollment_2019 %>% 
  pivot_wider(names_from='label', values_from=c('estimate'))

write.csv(
  public_enrollment_2019_wide,
  '~/epimodels/notebooks/AustinGranularModel/Schools/data/2019_AISD_Public_School_Enrollment_Estimates_by_ZCTA.csv'
)

write.csv(
  private_enrollment_2019_wide,
  '~/epimodels/notebooks/AustinGranularModel/Schools/data/2019_AISD_Private_School_Enrollment_Estimates_by_ZCTA.csv'
)