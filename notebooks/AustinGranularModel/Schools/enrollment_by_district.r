# school district populations
rm(list=ls())
library(tidyverse)
library(tidycensus)
library(sf)

API_KEY <- Sys.getenv("CENSUS_API_KEY")
census_api_key(API_KEY)

school_distr_pop = get_acs(
  geography='school district (unified)',
  variables=c('B09001_001', 'S0601_C01_001'),
  state='TX',
  year=2019,
  survey='acs5',
  show_call=TRUE,
  geometry=FALSE
)

school_distr_shape <- get_acs(
  geography='school district (unified)',
  variables=c('B09001_001'),
  state='TX',
  year=2019,
  survey='acs5',
  show_call=TRUE,
  geometry=TRUE
)

school_distr_pop <- school_distr_pop %>% 
  pivot_wider(names_from='variable', values_from=c('estimate', 'moe'))

# make estimate, moe names interpretable
names(school_distr_pop) <- c('GEOID', 'NAME', 'POP_17UNDER', 'TOTAL_POP', 'MOE_POP_17UNDER', 'MOE_TOTAL_POP')

# drop unnecessary columns
school_distr_shape <- school_distr_shape %>% select(-estimate, -moe)

# combine the shapefile and the population data
school_distr_data <- inner_join(school_distr_shape, school_distr_pop, by=c('GEOID', 'NAME'))

# load the enrollment data
enrollment <- read.csv('~/Downloads/Enrollment Report_Statewide_Districts_Grade_2020-2021.csv', skip=4)

# filter notes at bottom of file
# the notes read "Ranges (e.g. <10 and <20) indicate counts are not available (masked) to comply with 
# the Family Educational Rights and Privacy Act (FERPA). Masked numbers are typically small although larger 
# numbers may be masked to prevent imputation.)
enrollment <- enrollment %>% filter(GRADE != '')

# masked values will have a "<" sign in the enrollment value; these must be converted to numeric
# type before summation. we will just remove the "<" and convert to an integer, rather than attempt
# to impute the actual enrollment.
rm_less_than <- function(x) return(as.numeric(sub("<", "", x)))
enrollment$ENROLLMENT <- sapply(enrollment$ENROLLMENT, rm_less_than)

# get total district-wide enrollment
total_enrollment <- enrollment %>% 
  group_by(YEAR, REGION, COUNTY.NAME, DISTRICT, DISTRICT.NAME, CHARTER.STATUS) %>% 
  summarise(ENROLLMENT = sum(ENROLLMENT))

# convert the census data district names to upper case with acronyms
to_cisd <- function(x) return(sub('Consolidated Independent School District, Texas', 'CISD', x))
to_isd <- function(x) return(sub('Independent School District, Texas', 'ISD', x))
school_distr_data$NAME <- sapply(school_distr_data$NAME, to_cisd)
school_distr_data$NAME <- sapply(school_distr_data$NAME, to_isd)
school_distr_data$NAME <- sapply(school_distr_data$NAME, str_to_upper)

# join enrollment data with census-derived school data
final_data <- st_as_sf(left_join(total_enrollment, school_distr_data, by=c('DISTRICT.NAME'='NAME')))

# derive some additional data
final_data$PCT_ENROLLED_UNDER17 <- final_data$ENROLLMENT/final_data$POP_17UNDER
sensible_data <- final_data %>% filter(PCT_ENROLLED_UNDER17 <= 1)
austin_msa <- final_data %>% filter(COUNTY.NAME %in% c('TRAVIS COUNTY', 'WILLIAMSON COUNTY', 'HAYS COUNTY', 'CALDWELL COUNTY', 'BASTROP COUNTY'))
austin_districts <- unique(austin_msa$GEOID)

plot(austin_msa['PCT_ENROLLED_UNDER17'])

st_write(final_data, '~/COVID19/austin-spatial-prev/2020_2021_TX_School_Enrollment.shp', driver='ESRI Shapefile')

v19 <- load_variables(2019, dataset='acs5/subject')
school_enrollment <- v19 %>% filter(grepl('ENROLL', v19$concept)) %>% filter(grepl('Percent', label)) # %>% filter(!('graduate' %in% 'label'))

enrollment_shape <- get_acs(
  geography='school district (unified)',
  variables=school_enrollment$name,
  state='TX',
  year=2019,
  survey='acs5',
  show_call=TRUE,
  geometry=FALSE
)
austin_enrollment <- enrollment_shape %>% filter(GEOID %in% austin_districts)

enrollment_shape_named <- left_join(austin_enrollment, school_enrollment, by=c('variable'='name'))
enrollment_shape_named <- enrollment_shape_named %>% select(-c(moe, variable))

private_enrollment <- enrollment_shape_named %>% filter(grepl('private', label))
public_enrollment <- enrollment_shape_named %>% filter(grepl('public', label))

public_enrollment_wide <- public_enrollment %>% 
  pivot_wider(names_from='label', values_from=c('estimate'))
long_names <- names(public_enrollment_wide)

clean_columns <- c(
  "GEOID",
  "NAME",
  "concept",
  "3up_enr",
  "3up_prek",
  "3up_k12",
  "3up_k",
  "3up_1-4",
  "3up_5-8",
  "3up_9-12",
  "3up_ugr",
  "3up_gr",
  "ug_gr",
  "m_ug_gr",
  "f_ug_gr",
  "3-4y",
  "3-4y_enr",
  "5-9y",
  "5-9rs_enr",
  "10-14y",
  "10-14y_enr",
  "15-17",
  "15-17_enr",
  "18-19y",
  "18-19y_enr",
  "20-24y",
  "20-24y_enr",
  "25-34y",
  "25-34y_enr",
  "35up",
  "35up_enr",
  "18-24y",
  "18-24y_ug_gr",
  "18-24y_m",
  "18-24y_m_ug_gr",
  "18-24y_f",
  "18-24y_f_ug_gr"
)
names(public_enrollment_wide) <- clean_columns

private_enrollment_wide <- private_enrollment %>% 
  pivot_wider(names_from='label', values_from=c('estimate'))
long_names_private <- names(private_enrollment_wide)
names(private_enrollment_wide) <- clean_columns

final_austin <- st_as_sf(left_join(public_enrollment_wide, school_distr_shape, by=c('GEOID'='GEOID')))
final_austin_private <- st_as_sf(left_join(private_enrollment_wide, school_distr_shape, by=c('GEOID'='GEOID')))

st_write(final_austin, '~/COVID19/austin-spatial-prev/2019_AustinMSA_Census_Public_School_Enrollment_Estimates.shp', driver='ESRI Shapefile')

st_write(final_austin_private, '~/COVID19/austin-spatial-prev/2019_AustinMSA_Census_Private_School_Enrollment_Estimates.shp', driver='ESRI Shapefile')

column_decode <- tibble(
  clean_label=clean_columns,
  public_label=long_names,
  private_label=long_names_private
)

write.csv(column_decode, '~/COVID19/austin-spatial-prev/2019_AustinMSA_Census_Private_Public_School_Enrollment_Columns.csv')


school_enrollment_totals <- v19 %>% filter(grepl('ENROLL', v19$concept)) %>% filter(grepl('Total', label)) # %>% filter(!('graduate' %in% 'label'))
enrollment_totals_shape <- get_acs(
  geography='school district (unified)',
  variables=school_enrollment_totals$name,
  state='TX',
  year=2019,
  survey='acs5',
  show_call=TRUE,
  geometry=FALSE
)

austin_enrollment_totals <- enrollment_totals_shape %>% filter(GEOID %in% austin_districts)

enrollment_totals_shape_named <- left_join(austin_enrollment_totals, school_enrollment_totals, by=c('variable'='name'))
enrollment_totals_shape_named <- enrollment_totals_shape_named %>% select(-c(moe, variable))

library(stringr)

enrollment_totals_shape_named$join_name <- sapply(
  enrollment_totals_shape_named$label,
  sub, pattern='Total', replacement='__')

private_enrollment$join_name <- sapply(
  private_enrollment$label,
  sub, pattern='Percent in private school', replacement='__')

pct_total <- left_join(
  enrollment_totals_shape_named,
  private_enrollment,
  by=c('join_name', 'GEOID', 'NAME', 'concept'),
  suffix=c('_total', '_percent')
)

pct_total$estimate_fraction <- pct_total$estimate_percent / 100.0
pct_total$estimate_private <- pct_total$estimate_total * pct_total$estimate_fraction
pct_total$estimate_public <- pct_total$estimate_total - pct_total$estimate_private

percent_only <- pct_total %>% 
  select(GEOID, NAME, label_percent, estimate_percent) %>%
  mutate(measure_type='percent')
names(percent_only) <- c('GEOID', 'NAME', 'label', 'estimate', 'measure_type')

private_only <- pct_total %>% 
  select(GEOID, NAME, label_total, estimate_private) %>%
  mutate(measure_type='total_private')
names(private_only) <- c('GEOID', 'NAME', 'label', 'estimate', 'measure_type')

public_only <- pct_total %>% 
  select(GEOID, NAME, label_total, estimate_public) %>%
  mutate(measure_type='total_public')
names(public_only) <- c('GEOID', 'NAME', 'label', 'estimate', 'measure_type')

total_only <- pct_total %>% 
  select(GEOID, NAME, label_total, estimate_total) %>%
  mutate(measure_type='total_private_public')
names(total_only) <- c('GEOID', 'NAME', 'label', 'estimate', 'measure_type')

final_pct_total <- bind_rows(percent_only, total_only, public_only, private_only)
write.csv(final_pct_total, '~/COVID19/austin-spatial-prev/2019_AustinMSA_Census_Private_Public_School_Enrollment/2019_Austin_Percent_and_Total_Private_Enrollment.csv')

private_enrollment_total_wide <- private_enrollment_total %>% 
  pivot_wider(names_from='label', values_from=c('estimate'))
long_names_private_total <- names(private_enrollment_total_wide)
names(private_enrollment_total_wide) <- clean_columns
