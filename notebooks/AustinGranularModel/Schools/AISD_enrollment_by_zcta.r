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
school_enrollment_2018_percent <- v18 %>% filter(grepl('ENROLL', v18$concept)) %>% filter(grepl('Percent', label)) # %>% filter(!('graduate' %in% 'label'))
school_enrollment_2018_total <- v18 %>% filter(grepl('ENROLL', v18$concept)) %>% filter(grepl('Total', label)) # %>% filter(!('graduate' %in% 'label'))
school_enrollment_2018 <- bind_rows(school_enrollment_2018_percent, school_enrollment_2018_total)

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

write.csv(
  enrollment_2018_named,
  '~/epimodels/notebooks/AustinGranularModel/Schools/data/2018_AISD_Public_Private_School_Enrollment_Percent_and_Total_Estimates_by_ZCTA.csv'
  )

# repeat for 2019
# collect the enrollment variables and download enrollment data for AISD ZCTAs
v19 <- load_variables(2019, dataset='acs5/subject')
school_enrollment_2019_percent <- v19 %>% filter(grepl('ENROLL', v19$concept)) %>% filter(grepl('Percent', label)) # %>% filter(!('graduate' %in% 'label'))
school_enrollment_2019_total <- v19 %>% filter(grepl('ENROLL', v19$concept)) %>% filter(grepl('Total', label)) # %>% filter(!('graduate' %in% 'label'))
school_enrollment_2019 <- bind_rows(school_enrollment_2019_percent, school_enrollment_2019_total)

# geography hierarchies changed in 2019 so we have to retrieve data for all the ZCTAs...
enrollment_2019 <- get_acs(
  geography='zcta',
  variables=school_enrollment_2019$name,
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

write.csv(
  enrollment_2019_named,
  '~/epimodels/notebooks/AustinGranularModel/Schools/data/2019_AISD_Public_Private_School_Enrollment_Percent_and_Total_Estimates_by_ZCTA.csv'
)

### Combine data sets and extract select variables

enrollment_2018_named$year <- 2018
enrollment_2019_named$year <- 2019
enrollment_named <- bind_rows(enrollment_2018_named, enrollment_2019_named)

enrollment_vars <- c(
  "Estimate!!Total!!Population 5 to 9 years",                                                                                                                
  "Estimate!!Total!!Population 5 to 9 years!!5 to 9 year olds enrolled in school",                                                                            
  "Estimate!!Total!!Population 10 to 14 years",                                                                            
  "Estimate!!Total!!Population 10 to 14 years!!10 to 14 year olds enrolled in school",
  "Estimate!!Total!!Population 15 to 17",
  "Estimate!!Total!!Population 15 to 17!!15 to 17 year olds enrolled in school",
  "Estimate!!Percent in private school!!Population 5 to 9 years!!5 to 9 year olds enrolled in school",
  "Estimate!!Percent in private school!!Population 10 to 14 years!!10 to 14 year olds enrolled in school",
  "Estimate!!Percent in private school!!Population 15 to 17!!15 to 17 year olds enrolled in school"
)

enrollment_named_gradeschool <- enrollment_named %>%
  filter(label %in% enrollment_vars)

enrollment_gradeschool_wide <- enrollment_named_gradeschool %>% 
  pivot_wider(names_from='label', values_from=c('estimate'))

write.csv(enrollment_gradeschool_wide, 
          '~/epimodels/notebooks/AustinGranularModel/Schools/data/2018_2019_AISD_Select_Private_School_Enrollment_Estimates_by_ZCTA.csv')


