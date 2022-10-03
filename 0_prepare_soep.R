library(tidyverse)

# Format SOEP data
soep_filtered = read_csv("../../Data/SOEP/soep_gross_net_income_working_age_singles.csv") %>%
    mutate(
        wage = if_else(lfs_working == 1, monthly_wage_gross, 0),
        netwage = if_else(lfs_working == 1, monthly_wage_net, 0),
        part = (wage > 0),
        x0 = 1
    ) %>%
    filter(wage < 5000,!is.na(years_of_educ),!is.na(female),!is.na(single)) |>
    rename(
        x1 = years_of_educ,
        x2 = female,
        x3 = single,
        # x4 = age,
        y = wage
    ) |> 
    mutate(y0 = 0, w = 1) # No basic income in baseline data

soep_filtered |>
    write_csv("Data_processed/SOEP_with_income.csv")

soep_filtered |>
    select(x0, x1, x2, x3) |>
    write_csv("Data_processed/covariates_SOEP.csv")
