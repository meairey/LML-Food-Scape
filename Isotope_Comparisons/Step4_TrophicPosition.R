library(tidyverse)
library(broom)
#### Load in all isotope data then convert strings to numerics
load("Data/IsotopeComparisons/IsotopeData.RData")
full = data 



full_corrected = full %>%
  left_join(baseline_simmr %>% filter(GROUP == "B"), by = c("period" = "community")) %>%
  mutate(D13C_c = D13C - mean_d13C,
         D15N_c = (D15N - mean_d15N)/3.4 + 1)


# Filter your data first
dat_with_residuals <- full_corrected %>%
  filter(!is.na(D15N_c), !is.na(LENGTH)) %>%
  group_by(TAXON) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(D15N_c ~ log(LENGTH), data = .x)),
    augmented = map2(model, data, ~ augment(.x, data = .y))
  ) %>%
  select(TAXON, augmented) %>%
  unnest(augmented) %>%
  filter(TAXON %nin% c("BND","NRD"))

# Test: compare residuals between periods for each species
sig_results <- dat_with_residuals %>%
  filter(TAXON %nin% c("BND", "NRD", "LLS", "ST", "BB", "WS")) %>%
  group_by(TAXON) %>%
  summarize(
    p_value = wilcox.test(.resid ~ period)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )

# Find max δ15N value for each species
y_positions <- dat_with_residuals %>%
  group_by(TAXON) %>%
  summarize(
    max_y = max(D15N_c, na.rm = TRUE),
    .groups = "drop"
  )

# Combine with significance
sig_results <- sig_results %>%
  left_join(y_positions, by = "TAXON") %>%
  mutate(
    label_y = max_y + 0.15   # You can adjust the offset (0.2) as needed
  ) %>%
  left_join(suc, by = c("TAXON" = "CODE")) %>%
  left_join(legend_sci, by = c("TAXON" = "code"))
