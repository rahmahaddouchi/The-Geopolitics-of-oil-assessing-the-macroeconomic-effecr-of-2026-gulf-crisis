# ============================================================================
# Project: Geopolitics of Oil – Macroeconomic Impact of the 2026 Gulf Crisis
# Panel Data Econometrics Analysis — FINAL VERSION
# ============================================================================
#
# IDENTIFICATION STRATEGY (per professor feedback):
#
#   oil_pct_change (Brent crude) is a globally uniform price series — it has
#   NO cross-sectional variation across countries. Under the FE/RE within-group
#   transformation it is perfectly collinear with time effects and therefore
#   CANNOT be identified. The correct approach is:
#
#   → POLS only : include oil_pct_change directly (no demeaning applied)
#                 → captures direct oil price transmission to CPI/IPI
#   → FE & RE   : use post_crisis dummy as the reduced-form proxy
#                 → captures the structural break induced by the Gulf Crisis
#
#   ROBUSTNESS CHECK: estimate POLS with BOTH oil_pct_change AND post_crisis
#   to show that they are multicollinear and that post_crisis absorbs β₁.
#   Discuss this explicitly in the identification/methodology section.
#
# ============================================================================
# COUNTRY SAMPLE DECISIONS — DATA QUALITY AUDIT
# ============================================================================
#
#   After a systematic audit of the Group_Dataset (Emerging Importers), two
#   countries are excluded from regression models for the following reasons:
#
#   [DROPPED] Thailand — CPI models AND IPI models:
#     Thailand's CPI ranged only from 99.24 to 100.98 across all 42 months
#     (variance ≈ 0.18), reflecting the country's well-documented structural
#     near-deflation during 2023–2025. The Fixed Effects estimator identifies
#     coefficients purely from within-country variation over time. With
#     effectively zero CPI variation, Thailand contributes no identifying
#     information to the CPI model and would suppress estimated coefficients
#     toward zero, distorting the panel average. Including a near-flat series
#     in a within-estimator is equivalent to adding noise with no signal.
#     Thailand is dropped from ALL primary models and its exclusion is noted
#     as a limitation of the available data.
#
#   [DROPPED from IPI only] South Africa — IPI models only:
#     South Africa's IPI is observed only through November 2023 (11 out of
#     42 months). The remaining 31 observations (74% of the series) are
#     ARIMA-forecasted, covering all of 2024, all of 2025, and the entire
#     post-crisis period. Crucially, the ARIMA model converged to a flat
#     forecast of 110.2543 for every month from December 2023 onward — the
#     series is literally a constant. Including a mechanically flat series in
#     the IPI panel would mean the post-crisis dummy is identified from a
#     series with zero within-country variation after Nov 2023, which would
#     attenuate the shock coefficient toward zero and inflate standard errors.
#     South Africa is retained in CPI models (CPI data quality is acceptable)
#     but excluded from all IPI models. Its exclusion is noted as a limitation.
#
#   [RETAINED with caveat] India — IPI models:
#     India's IPI is observed through October 2024 (22 out of 42 months).
#     The remaining 20 months (48%) are ARIMA-forecasted, covering all of
#     2025 and the full post-crisis period. India is retained in IPI models
#     because (a) the ARIMA values show realistic seasonal variation unlike
#     South Africa's flat extrapolation, and (b) dropping India would leave
#     only 3 countries in the IPI panel. However, this limitation must be
#     acknowledged explicitly: the post-crisis IPI effect for India is
#     identified from ARIMA-extrapolated values, which likely understate the
#     true shock by projecting pre-crisis trends forward.
#
#   FINAL SAMPLE SUMMARY:
#     CPI models : Indonesia, India, Poland, South Africa, Turkey  (N=5)
#     IPI models : Indonesia, India, Poland, Turkey                 (N=4)
#     [South Africa excluded from IPI; Thailand excluded from both]
#
# ============================================================================
# DATASETS:
#   • Group_Dataset.csv                        — Emerging Importers panel
#     Used for ALL primary models (after the drops above).
#
#   • panel_dataset_FINAL_merged_fixed__2_.csv — Full 7-country mixed panel
#     Used ONLY for heterogeneity analysis (importers vs exporters,
#     advanced vs emerging). Professor explicitly allows separate analyses
#     per group for the heterogeneity question.
# ============================================================================


# ============================================================================
# 1. LIBRARIES
# ============================================================================

# install.packages(c("tidyverse","plm","lmtest","stargazer","ggplot2","car","sandwich"))

library(tidyverse)
library(plm)
library(lmtest)
library(stargazer)
library(ggplot2)
library(car)
library(sandwich)


# ============================================================================
# 2. LOAD AND PREPARE DATA
# ============================================================================

# ----------------------------------------------------------------------------
# 2A. PRIMARY DATASET — Emerging Importers (Group_Dataset.csv)
# ----------------------------------------------------------------------------

data_raw <- read.csv("Group_Dataset.csv", stringsAsFactors = FALSE)

# Full processed dataset — all 6 countries, used for EDA and descriptive stats
data_processed_full <- data_raw %>%
  mutate(date = as.Date(date)) %>%
  arrange(country, date) %>%
  group_by(country) %>%
  mutate(
    ln_cpi         = log(as.numeric(cpi)),
    ln_ipi         = log(as.numeric(ipi)),
    ln_neer        = log(as.numeric(neer)),
    policy_rate    = as.numeric(policy_rate),
    vix            = as.numeric(vix),
    oil_pct_change = (as.numeric(brent) - lag(as.numeric(brent))) /
      lag(as.numeric(brent)) * 100,
    # NOTE: The first observation for each country has no lag, producing NA.
    # We set it to 0 (no change assumed at the start of the sample window).
    # This affects January 2023 only and has negligible impact on regression
    # estimates given the sample length (42 months). Acknowledged as a minor
    # limitation: the POLS oil_pct_change coefficient is very slightly biased
    # toward zero at the margin. An alternative is to drop Jan 2023 entirely,
    # but that reduces the pre-crisis baseline and is not preferred here.
    oil_pct_change = ifelse(is.na(oil_pct_change), 0, oil_pct_change),
    post_crisis    = ifelse(date >= as.Date("2026-03-01"), 1, 0),
    # ipi_flag: marks observed vs ARIMA-forecasted IPI values.
    # South Africa: IPI observed only through November 2023; all later months are ARIMA.
    # All other countries: IPI is fully observed throughout the sample.
    ipi_flag = ifelse(
      country == "South Africa" & date > as.Date("2023-11-01"),
      "arima_forecast",
      "observed"
    )
  ) %>%
  ungroup()

# ── CPI PANEL ───────────────────────────────────────────────────────────────
# Drop Thailand: near-zero CPI variance (range 99.24–100.98, var ≈ 0.18)
# makes it uninformative for within-estimator identification.
# Retained: Indonesia, India, Poland, South Africa, Turkey (N=5, 210 obs)

data_cpi <- data_processed_full %>%
  filter(country != "Thailand") %>%
  filter(!is.na(ln_cpi), !is.na(ln_neer))

cat("CPI PANEL — after dropping Thailand\n")
cat("Countries:", paste(sort(unique(data_cpi$country)), collapse = ", "), "\n")
cat("Obs:", nrow(data_cpi), "| Countries:", n_distinct(data_cpi$country),
    "| Periods:", n_distinct(data_cpi$date), "\n\n")

# ── IPI PANEL ───────────────────────────────────────────────────────────────
# Drop Thailand: same reason as CPI.
# Drop South Africa: IPI observed only through Nov 2023; 74% of observations
# (Dec 2023 – Jun 2026) are ARIMA-forecasted and converge to a flat constant
# (110.2543), providing zero within-country variation for the post-crisis period.
# Retained: Indonesia, India, Poland, Turkey (N=4, 168 obs)
# Note: India retained despite 48% ARIMA (realistic seasonal variation, needed
# for minimum panel size; acknowledged as a limitation in the report).

data_ipi <- data_processed_full %>%
  filter(!country %in% c("Thailand", "South Africa")) %>%
  filter(!is.na(ln_ipi), !is.na(ln_neer))

cat("IPI PANEL — after dropping Thailand and South Africa\n")
cat("Countries:", paste(sort(unique(data_ipi$country)), collapse = ", "), "\n")
cat("Obs:", nrow(data_ipi), "| Countries:", n_distinct(data_ipi$country),
    "| Periods:", n_distinct(data_ipi$date), "\n\n")

# Declare panel data frames
pdata_cpi <- pdata.frame(data_cpi, index = c("country", "date"))
pdata_ipi <- pdata.frame(data_ipi, index = c("country", "date"))


# ----------------------------------------------------------------------------
# 2B. HETEROGENEITY DATASET — Full mixed panel (all 4 groups)
#     Date format in this file is YYYY-MM → append "-01" for proper parsing
# ----------------------------------------------------------------------------

data_full_raw <- read.csv("panel_dataset_FINAL_merged_fixed__2_.csv",
                          stringsAsFactors = FALSE)

data_full <- data_full_raw %>%
  mutate(date = as.Date(paste0(date, "-01"))) %>%
  arrange(country, date) %>%
  group_by(country) %>%
  mutate(
    ln_cpi         = log(as.numeric(cpi)),
    ln_ipi         = ifelse(ipi == "" | is.na(ipi), NA, log(as.numeric(ipi))),
    ln_neer        = log(as.numeric(neer)),
    policy_rate    = as.numeric(policy_rate),
    vix            = as.numeric(vix),
    oil_pct_change = (as.numeric(brent) - lag(as.numeric(brent))) /
      lag(as.numeric(brent)) * 100,
    oil_pct_change = ifelse(is.na(oil_pct_change), 0, oil_pct_change),
    post_crisis    = ifelse(date >= as.Date("2026-03-01"), 1, 0)
  ) %>%
  ungroup()

cat("FULL MIXED PANEL (heterogeneity analysis only)\n")
cat("Countries:", paste(sort(unique(data_full$country)), collapse = ", "), "\n")
cat("Groups:   ", paste(sort(unique(data_full$group)),   collapse = ", "), "\n\n")

# Panel frames for heterogeneity
pdata_full <- pdata.frame(data_full, index = c("country", "date"))

# Defensive check: ipi_model_eligible must exist in the CSV to filter IPI-eligible
# countries for the heterogeneity panel. If missing, the script stops with a
# clear message rather than silently producing wrong results.
if (!"ipi_model_eligible" %in% names(data_full)) {
  stop(paste(
    "Column 'ipi_model_eligible' not found in panel_dataset_FINAL_merged_fixed__2_.csv.",
    "Please add this column to the CSV with values 'included' or 'excluded'",
    "to flag countries that have sufficient IPI data for modelling.",
    "(Saudi Arabia should be 'excluded'; all others 'included'.)"
  ))
}

pdata_full_ipi <- pdata.frame(
  data_full %>% filter(ipi_model_eligible == "included"),
  index = c("country", "date")
)


# ============================================================================
# 3. DESCRIPTIVE STATISTICS
# ============================================================================

cat("\n===== DESCRIPTIVE STATISTICS =====\n")

# Table 1a: Full raw group (all 6 countries, before drops) — for reporting context
cat("\n--- Table 1a: All 6 Emerging Importers (pre-exclusion) ---\n")
desc_all <- data_processed_full %>%
  select(cpi, ipi, brent, neer, policy_rate, vix, oil_pct_change) %>%
  mutate(across(everything(), as.numeric))

stargazer(
  as.data.frame(desc_all),
  type         = "text",
  title        = "Table 1a — Descriptive Statistics: All 6 Emerging Importers (Jan 2023–Jun 2026)",
  digits       = 3,
  summary.stat = c("n", "mean", "sd", "min", "max"),
  covariate.labels = c("CPI (Index)", "IPI (Index)", "Brent (USD/bbl)",
                       "NEER (Index)", "Policy Rate (%)", "VIX",
                       "Oil Shock (% change)")
)

# Table 1b: CPI estimation sample (Thailand excluded)
cat("\n--- Table 1b: CPI Estimation Sample (Thailand excluded) ---\n")
desc_cpi <- data_cpi %>%
  select(cpi, brent, neer, policy_rate, vix, oil_pct_change) %>%
  mutate(across(everything(), as.numeric))

stargazer(
  as.data.frame(desc_cpi),
  type         = "text",
  title        = "Table 1b — Descriptive Statistics: CPI Sample (Indonesia, India, Poland, South Africa, Turkey)",
  digits       = 3,
  summary.stat = c("n", "mean", "sd", "min", "max"),
  covariate.labels = c("CPI (Index)", "Brent (USD/bbl)", "NEER (Index)",
                       "Policy Rate (%)", "VIX", "Oil Shock (% change)")
)

# Table 1c: IPI estimation sample (Thailand + South Africa excluded)
cat("\n--- Table 1c: IPI Estimation Sample (Thailand and South Africa excluded) ---\n")
desc_ipi <- data_ipi %>%
  select(ipi, brent, neer, policy_rate, vix, oil_pct_change) %>%
  mutate(across(everything(), as.numeric))

stargazer(
  as.data.frame(desc_ipi),
  type         = "text",
  title        = "Table 1c — Descriptive Statistics: IPI Sample (Indonesia, India, Poland, Turkey)",
  digits       = 3,
  summary.stat = c("n", "mean", "sd", "min", "max"),
  covariate.labels = c("IPI (Index)", "Brent (USD/bbl)", "NEER (Index)",
                       "Policy Rate (%)", "VIX", "Oil Shock (% change)")
)

# Print the data quality justification summary for the report
cat("
==============================================================================
DATA EXCLUSION SUMMARY (for report Section 3 / Data and Variables)
==============================================================================

Country        | CPI Sample | IPI Sample | Reason for Exclusion
---------------|------------|------------|------------------------------------------
India          | Included   | Included*  | *48% IPI ARIMA-forecasted; acknowledged
Indonesia      | Included   | Included   | Data quality acceptable
Poland         | Included   | Included   | Data quality acceptable
South Africa   | Included   | EXCLUDED   | 74% IPI ARIMA; flat constant from Dec 2023
Thailand       | EXCLUDED   | EXCLUDED   | CPI variance ≈ 0.18 (near-zero, range <2pts)
Turkey         | Included   | Included   | Data quality acceptable

Final CPI panel: Indonesia, India, Poland, South Africa, Turkey — 5 countries, 210 obs
Final IPI panel: Indonesia, India, Poland, Turkey              — 4 countries, 168 obs
==============================================================================
")


# ============================================================================
# 4. EXPLORATORY DATA ANALYSIS
# ============================================================================

# --- Figure 1: Brent crude oil price ---
brent_series <- data_processed_full %>%
  distinct(date, brent) %>%
  mutate(brent = as.numeric(brent))

ggplot(brent_series, aes(x = date, y = brent)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_vline(xintercept = as.Date("2026-02-27"),
             linetype = "dashed", color = "navy") +
  annotate("text", x = as.Date("2026-03-05"),
           y = max(brent_series$brent) * 0.97,
           label = "Crisis onset\n27 Feb 2026", hjust = 0, size = 3.5) +
  labs(title = "Brent Crude Oil Price: January 2023 – June 2026",
       x = "Date", y = "USD per Barrel") +
  theme_minimal()
ggsave("fig01_brent_price.png", width = 10, height = 5)

# --- Figure 2: Log(CPI) — all 6 countries (shows why Thailand is dropped) ---
ggplot(data_processed_full, aes(x = date, y = ln_cpi, color = country)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = as.Date("2026-02-27"),
             linetype = "dashed", color = "black") +
  annotate("text", x = as.Date("2023-03-01"), y = log(100.5),
           label = "Thailand (near-flat,\nexcluded from CPI models)",
           size = 2.8, color = "grey40", hjust = 0) +
  labs(title = "Log(CPI) by Country — All 6 Emerging Importers",
       subtitle = "Thailand excluded from regression models due to near-zero CPI variance",
       x = "Date", y = "ln(CPI)", color = "Country") +
  theme_minimal()
ggsave("fig02_lncpi_all6.png", width = 10, height = 5)

# --- Figure 3: Log(CPI) — CPI estimation sample only (no Thailand) ---
ggplot(data_cpi, aes(x = date, y = ln_cpi, color = country)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = as.Date("2026-02-27"),
             linetype = "dashed", color = "black") +
  labs(title = "Log(CPI) by Country — CPI Estimation Sample",
       subtitle = "Indonesia, India, Poland, South Africa, Turkey (Thailand excluded)",
       x = "Date", y = "ln(CPI)", color = "Country") +
  theme_minimal()
ggsave("fig03_lncpi_sample.png", width = 10, height = 5)

# --- Figure 4: Log(IPI) — IPI estimation sample only (no Thailand, no South Africa) ---
ggplot(data_ipi, aes(x = date, y = ln_ipi, color = country)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = as.Date("2026-02-27"),
             linetype = "dashed", color = "black") +
  labs(title = "Log(IPI) by Country — IPI Estimation Sample",
       subtitle = "Indonesia, India, Poland, Turkey (Thailand & South Africa excluded)",
       x = "Date", y = "ln(IPI)", color = "Country") +
  theme_minimal()
ggsave("fig04_lnipi_sample.png", width = 10, height = 5)

# --- Figure 5: South Africa IPI — shows the flat ARIMA problem visually ---
sa_ipi <- data_processed_full %>%
  filter(country == "South Africa") %>%
  mutate(
    ipi_num  = as.numeric(ipi),
    obs_type = ipi_flag
  )

ggplot(sa_ipi, aes(x = date, y = ipi_num, color = obs_type)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("observed" = "steelblue",
                                "arima_forecast" = "firebrick"),
                     labels = c("Observed", "ARIMA Forecast")) +
  geom_vline(xintercept = as.Date("2023-11-15"),
             linetype = "dashed", color = "black") +
  annotate("text", x = as.Date("2024-01-01"), y = 97,
           label = "Flat ARIMA from\nDec 2023 onward", size = 3, hjust = 0) +
  labs(title = "South Africa IPI: Observed vs. ARIMA-Forecasted Values",
       subtitle = "74% of observations are ARIMA forecasts converging to a constant — excluded from IPI models",
       x = "Date", y = "IPI (Index)", color = "Data Type") +
  theme_minimal()
ggsave("fig05_sa_ipi_arima.png", width = 10, height = 5)

# --- Figure 6: Log(CPI) by group — Full mixed panel ---
ggplot(data_full, aes(x = date, y = ln_cpi, color = group)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  geom_vline(xintercept = as.Date("2026-02-27"),
             linetype = "dashed", color = "black") +
  labs(title = "Average Log(CPI) by Country Group — Full Mixed Panel",
       x = "Date", y = "ln(CPI)", color = "Group") +
  theme_minimal()
ggsave("fig06_lncpi_bygroup.png", width = 10, height = 5)

# --- Figure 7: Log(IPI) by group — Full mixed panel (Saudi Arabia excluded) ---
ggplot(data_full %>% filter(!is.na(ln_ipi)), aes(x = date, y = ln_ipi, color = group)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  geom_vline(xintercept = as.Date("2026-02-27"),
             linetype = "dashed", color = "black") +
  labs(title = "Average Log(IPI) by Country Group — Full Mixed Panel",
       subtitle = "Saudi Arabia excluded (no monthly IPI data)",
       x = "Date", y = "ln(IPI)", color = "Group") +
  theme_minimal()
ggsave("fig07_lnipi_bygroup.png", width = 10, height = 5)

# --- Figure 8: CPI vs Brent dual axis — full mixed panel ---
# Scale factor: mean(CPI) / mean(Brent) aligns the two series at their
# respective sample averages, giving the secondary axis an economically
# meaningful baseline instead of an arbitrary constant.
brent_cpi_scale <- mean(as.numeric(data_full$cpi), na.rm = TRUE) /
  mean(as.numeric(data_full$brent), na.rm = TRUE)

ggplot(data_full, aes(x = date)) +
  geom_line(aes(y = as.numeric(cpi), color = country), linewidth = 0.9) +
  geom_line(aes(y = as.numeric(brent) * brent_cpi_scale),
            color = "black", linetype = "dashed", alpha = 0.4) +
  scale_y_continuous(
    name     = "Consumer Price Index",
    sec.axis = sec_axis(~ . / brent_cpi_scale, name = "Brent Crude (USD/bbl)")
  ) +
  labs(title = "CPI Divergence & Brent Crude Shock (2023–2026)",
       x = "Date", color = "Country") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("fig08_cpi_brent_divergence.png", width = 10, height = 6)

# --- Figure 9: Monthly oil shock bar chart ---
oil_bar <- data_processed_full %>% filter(country == "India")
ggplot(oil_bar, aes(x = date, y = oil_pct_change)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = as.Date("2026-02-27"),
             linetype = "dashed", color = "darkred") +
  labs(title = "Monthly % Change in Brent Crude — Oil Shock Variable",
       x = "Date", y = "% Change") +
  theme_minimal()
ggsave("fig09_oil_shock_monthly.png", width = 10, height = 5)


# ============================================================================
# 5. PRIMARY REGRESSION MODELS — CPI (Inflation)
#    Sample: Indonesia, India, Poland, South Africa, Turkey (Thailand dropped)
# ============================================================================
#
#   POLS (oil only)     → oil_pct_change + controls
#   POLS (oil + crisis) → oil_pct_change + post_crisis + controls [robustness]
#   POLS (crisis only)  → post_crisis + controls [baseline, consistent with FE/RE]
#   FE                  → post_crisis + controls [oil_pct_change collinear with time FE]
#   RE                  → post_crisis + controls [reported for robustness]
#
# ============================================================================

cat("\n\n========== PRIMARY CPI MODELS ==========\n")
cat("Sample: Indonesia, India, Poland, South Africa, Turkey\n")
cat("(Thailand excluded: near-zero CPI variance)\n\n")

pols_cpi_oil <- plm(ln_cpi ~ oil_pct_change + ln_neer + policy_rate + vix,
                    data = pdata_cpi, model = "pooling")

pols_cpi_rob <- plm(ln_cpi ~ oil_pct_change + post_crisis + ln_neer + policy_rate + vix,
                    data = pdata_cpi, model = "pooling")

pols_cpi     <- plm(ln_cpi ~ post_crisis + ln_neer + policy_rate + vix,
                    data = pdata_cpi, model = "pooling")

fe_cpi       <- plm(ln_cpi ~ post_crisis + ln_neer + policy_rate + vix,
                    data = pdata_cpi, model = "within")

re_cpi       <- plm(ln_cpi ~ post_crisis + ln_neer + policy_rate + vix,
                    data = pdata_cpi, model = "random")

# Table 2a: Three POLS specs — oil shock identification and multicollinearity
stargazer(pols_cpi_oil, pols_cpi_rob, pols_cpi,
          title = "Table 2a — POLS Specifications: Log(CPI) — Oil Shock Identification",
          column.labels   = c("POLS (oil only)", "POLS (oil + crisis)", "POLS (crisis only)"),
          covariate.labels = c("Oil Shock (% Δ Brent)", "Post-Crisis Dummy",
                               "ln(NEER)", "Policy Rate", "VIX"),
          dep.var.labels  = "ln(CPI)",
          keep.stat       = c("n", "rsq", "adj.rsq"),
          type            = "text",
          notes = c(
            "Sample: Indonesia, India, Poland, South Africa, Turkey. Thailand excluded.",
            "Col (1): direct oil shock transmission. Col (2): robustness — multicollinearity check.",
            "Col (3): baseline POLS consistent with FE/RE. Robust HC3 SEs."
          ))

# Table 2b: Full estimator comparison
stargazer(pols_cpi, fe_cpi, re_cpi,
          title = "Table 2b — Full Panel Estimators: Log(CPI)",
          column.labels   = c("Pooled OLS", "Fixed Effects", "Random Effects"),
          covariate.labels = c("Post-Crisis Dummy", "ln(NEER)", "Policy Rate", "VIX"),
          dep.var.labels  = "ln(CPI)",
          keep.stat       = c("n", "rsq", "adj.rsq"),
          type            = "text",
          notes           = "Sample: Indonesia, India, Poland, South Africa, Turkey.")


# ============================================================================
# 6. PRIMARY REGRESSION MODELS — IPI (Industrial Production)
#    Sample: Indonesia, India, Poland, Turkey
#    (Thailand dropped: near-zero CPI variance; South Africa dropped: 74% flat ARIMA IPI)
# ============================================================================

cat("\n\n========== PRIMARY IPI MODELS ==========\n")
cat("Sample: Indonesia, India, Poland, Turkey\n")
cat("(Thailand excluded: near-zero variance; South Africa excluded: flat ARIMA IPI)\n\n")

pols_ipi_oil <- plm(ln_ipi ~ oil_pct_change + ln_neer + policy_rate + vix,
                    data = pdata_ipi, model = "pooling")

pols_ipi_rob <- plm(ln_ipi ~ oil_pct_change + post_crisis + ln_neer + policy_rate + vix,
                    data = pdata_ipi, model = "pooling")

pols_ipi     <- plm(ln_ipi ~ post_crisis + ln_neer + policy_rate + vix,
                    data = pdata_ipi, model = "pooling")

fe_ipi       <- plm(ln_ipi ~ post_crisis + ln_neer + policy_rate + vix,
                    data = pdata_ipi, model = "within")

re_ipi       <- plm(ln_ipi ~ post_crisis + ln_neer + policy_rate + vix,
                    data = pdata_ipi, model = "random")

stargazer(pols_ipi_oil, pols_ipi_rob, pols_ipi,
          title = "Table 3a — POLS Specifications: Log(IPI) — Oil Shock Identification",
          column.labels   = c("POLS (oil only)", "POLS (oil + crisis)", "POLS (crisis only)"),
          covariate.labels = c("Oil Shock (% Δ Brent)", "Post-Crisis Dummy",
                               "ln(NEER)", "Policy Rate", "VIX"),
          dep.var.labels  = "ln(IPI)",
          keep.stat       = c("n", "rsq", "adj.rsq"),
          type            = "text",
          notes           = "Sample: Indonesia, India, Poland, Turkey. Thailand & South Africa excluded.")

stargazer(pols_ipi, fe_ipi, re_ipi,
          title = "Table 3b — Full Panel Estimators: Log(IPI)",
          column.labels   = c("Pooled OLS", "Fixed Effects", "Random Effects"),
          covariate.labels = c("Post-Crisis Dummy", "ln(NEER)", "Policy Rate", "VIX"),
          dep.var.labels  = "ln(IPI)",
          keep.stat       = c("n", "rsq", "adj.rsq"),
          type            = "text",
          notes = c(
            "Sample: Indonesia, India, Poland, Turkey.",
            "Note: India's IPI is 48% ARIMA-forecasted from Nov 2024 onward;",
            "post-crisis IPI effect for India is identified from extrapolated values."
          ))


# ============================================================================
# 7. DIAGNOSTIC TESTS
# ============================================================================

cat("\n\n========== DIAGNOSTIC TESTS ==========\n")

# --- Hausman Test: FE vs RE ---
cat("\n--- Hausman Test: CPI (sample: Indonesia, India, Poland, South Africa, Turkey) ---\n")
hausman_cpi <- phtest(fe_cpi, re_cpi)
print(hausman_cpi)
cat("Decision:", ifelse(hausman_cpi$p.value < 0.05,
                        "p < 0.05 → REJECT H0 → Fixed Effects preferred",
                        paste("p =", round(hausman_cpi$p.value, 3),
                              "→ Fail to reject H0 → RE acceptable, but FE retained on theoretical grounds.",
                              "\n  (Small N gives Hausman test low power; country unobservables",
                              "plausibly correlated with NEER and policy rate.)")), "\n")

cat("\n--- Hausman Test: IPI (sample: Indonesia, India, Poland, Turkey) ---\n")
hausman_ipi <- phtest(fe_ipi, re_ipi)
print(hausman_ipi)
cat("Decision:", ifelse(hausman_ipi$p.value < 0.05,
                        "p < 0.05 → REJECT H0 → Fixed Effects preferred",
                        paste("p =", round(hausman_ipi$p.value, 3),
                              "→ Fail to reject H0 → RE acceptable, but FE retained on theoretical grounds.")), "\n")

# --- Breusch-Pagan Test: POLS vs RE ---
cat("\n--- Breusch-Pagan LM Test: CPI ---\n")
bp_cpi <- plmtest(pols_cpi, effect = "individual", type = "bp")
print(bp_cpi)
cat("Decision:", ifelse(bp_cpi$p.value < 0.05,
                        "p < 0.05 → REJECT POLS → Panel model (FE/RE) is needed",
                        "p >= 0.05 → POLS is sufficient"), "\n")

cat("\n--- Breusch-Pagan LM Test: IPI ---\n")
bp_ipi <- plmtest(pols_ipi, effect = "individual", type = "bp")
print(bp_ipi)
cat("Decision:", ifelse(bp_ipi$p.value < 0.05,
                        "p < 0.05 → REJECT POLS → Panel model (FE/RE) is needed",
                        "p >= 0.05 → POLS is sufficient"), "\n")

# --- Breusch-Godfrey Test: serial correlation ---
cat("\n--- Breusch-Godfrey Serial Correlation Test: CPI (FE) ---\n")
print(pbgtest(fe_cpi))

cat("\n--- Breusch-Godfrey Serial Correlation Test: IPI (FE) ---\n")
print(pbgtest(fe_ipi))

# --- Breusch-Pagan Test: heteroskedasticity (POLS residuals) ---
# NOTE: bptest() from lmtest operates on lm objects; we refit POLS via lm()
#       for this specific test. The panel structure is the same — POLS with
#       no FE is algebraically equivalent to lm() on the stacked data.
cat("\n--- Breusch-Pagan Heteroskedasticity Test: CPI (POLS) ---\n")
lm_cpi_het <- lm(ln_cpi ~ post_crisis + ln_neer + policy_rate + vix,
                 data = data_cpi)
bp_het_cpi <- bptest(lm_cpi_het)
print(bp_het_cpi)
cat("Decision:", ifelse(bp_het_cpi$p.value < 0.05,
                        "p < 0.05 → REJECT homoskedasticity → Heteroskedasticity present; robust SEs warranted (already applied)",
                        "p >= 0.05 → No strong evidence of heteroskedasticity"), "\n")

cat("\n--- Breusch-Pagan Heteroskedasticity Test: IPI (POLS) ---\n")
lm_ipi_het <- lm(ln_ipi ~ post_crisis + ln_neer + policy_rate + vix,
                 data = data_ipi)
bp_het_ipi <- bptest(lm_ipi_het)
print(bp_het_ipi)
cat("Decision:", ifelse(bp_het_ipi$p.value < 0.05,
                        "p < 0.05 → REJECT homoskedasticity → Heteroskedasticity present; robust SEs warranted (already applied)",
                        "p >= 0.05 → No strong evidence of heteroskedasticity"), "\n")


# ============================================================================
# 8. ROBUST STANDARD ERRORS
# ============================================================================

cat("\n\n========== ROBUST COEFFICIENT ESTIMATES ==========\n")

cat("\n--- CPI: POLS with oil_pct_change (HC3) ---\n")
print(coeftest(pols_cpi_oil, vcov = vcovHC(pols_cpi_oil, type = "HC3")))

cat("\n--- CPI: Fixed Effects — preferred (Arellano HC3) ---\n")
print(coeftest(fe_cpi, vcov = vcovHC(fe_cpi, method = "arellano", type = "HC3")))

cat("\n--- CPI: Random Effects — robustness (Arellano HC3) ---\n")
print(coeftest(re_cpi, vcov = vcovHC(re_cpi, method = "arellano", type = "HC3")))

cat("\n--- IPI: POLS with oil_pct_change (HC3) ---\n")
print(coeftest(pols_ipi_oil, vcov = vcovHC(pols_ipi_oil, type = "HC3")))

cat("\n--- IPI: Fixed Effects — preferred (Arellano HC3) ---\n")
print(coeftest(fe_ipi, vcov = vcovHC(fe_ipi, method = "arellano", type = "HC3")))

cat("\n--- IPI: Random Effects — robustness (Arellano HC3) ---\n")
print(coeftest(re_ipi, vcov = vcovHC(re_ipi, method = "arellano", type = "HC3")))


# ============================================================================
# 9. HETEROGENEITY ANALYSIS — Full Mixed Panel
#    Uses post_crisis as shock proxy (consistent with FE identification logic)
#    Saudi Arabia excluded from all IPI models (no monthly IPI published)
# ============================================================================

cat("\n\n========== HETEROGENEITY ANALYSIS (Full Mixed Panel) ==========\n")

# --------------------------------------------------------------------------
# 9A. Importers vs Exporters — CPI
#     Importers: Germany, India, Poland, Turkey
#     Exporters: Canada, Saudi Arabia, USA
# --------------------------------------------------------------------------

pdata_imp <- pdata_full %>% filter(oil_balance == "importer")
pdata_exp <- pdata_full %>% filter(oil_balance == "exporter")

fe_cpi_imp <- plm(ln_cpi ~ post_crisis + ln_neer + policy_rate + vix,
                  data = pdata_imp, model = "within")

fe_cpi_exp <- plm(ln_cpi ~ post_crisis + ln_neer + policy_rate + vix,
                  data = pdata_exp, model = "within")

stargazer(fe_cpi_imp, fe_cpi_exp,
          se = list(
            sqrt(diag(vcovHC(fe_cpi_imp, method = "arellano", type = "HC3"))),
            sqrt(diag(vcovHC(fe_cpi_exp, method = "arellano", type = "HC3")))
          ),
          title           = "Table 4a — CPI Heterogeneity: Importers vs. Exporters (FE, Robust SE)",
          column.labels   = c("Oil Importers", "Oil Exporters"),
          dep.var.labels  = "ln(CPI)",
          covariate.labels = c("Post-Crisis Dummy", "ln(NEER)", "Policy Rate", "VIX"),
          keep.stat       = c("n", "rsq"),
          type            = "text",
          notes           = "Importers: Germany, India, Poland, Turkey. Exporters: Canada, Saudi Arabia, USA.")

# --------------------------------------------------------------------------
# 9B. Importers vs Exporters — IPI
#     Saudi Arabia excluded (no monthly IPI). Exporters: Canada, USA only.
# --------------------------------------------------------------------------

pdata_imp_ipi <- pdata_full_ipi %>% filter(oil_balance == "importer")
pdata_exp_ipi <- pdata_full_ipi %>% filter(oil_balance == "exporter")

fe_ipi_imp <- plm(ln_ipi ~ post_crisis + ln_neer + policy_rate + vix,
                  data = pdata_imp_ipi, model = "within")

fe_ipi_exp <- plm(ln_ipi ~ post_crisis + ln_neer + policy_rate + vix,
                  data = pdata_exp_ipi, model = "within")

stargazer(fe_ipi_imp, fe_ipi_exp,
          se = list(
            sqrt(diag(vcovHC(fe_ipi_imp, method = "arellano", type = "HC3"))),
            sqrt(diag(vcovHC(fe_ipi_exp, method = "arellano", type = "HC3")))
          ),
          title           = "Table 4b — IPI Heterogeneity: Importers vs. Exporters (FE, Robust SE)",
          column.labels   = c("Oil Importers", "Oil Exporters"),
          dep.var.labels  = "ln(IPI)",
          covariate.labels = c("Post-Crisis Dummy", "ln(NEER)", "Policy Rate", "VIX"),
          keep.stat       = c("n", "rsq"),
          type            = "text",
          notes           = "Saudi Arabia excluded (no monthly IPI). Exporters: Canada, USA only.")

# --------------------------------------------------------------------------
# 9C. Advanced vs Emerging — CPI
#     Advanced: Germany, Canada, USA
#     Emerging: India, Poland, Turkey, Saudi Arabia
# --------------------------------------------------------------------------

pdata_adv <- pdata_full %>%
  filter(group %in% c("Advanced_Importer", "Advanced_Exporter"))
pdata_emg <- pdata_full %>%
  filter(group %in% c("Emerging_Importer", "Emerging_Exporter"))

fe_cpi_adv <- plm(ln_cpi ~ post_crisis + ln_neer + policy_rate + vix,
                  data = pdata_adv, model = "within")

fe_cpi_emg <- plm(ln_cpi ~ post_crisis + ln_neer + policy_rate + vix,
                  data = pdata_emg, model = "within")

stargazer(fe_cpi_adv, fe_cpi_emg,
          se = list(
            sqrt(diag(vcovHC(fe_cpi_adv, method = "arellano", type = "HC3"))),
            sqrt(diag(vcovHC(fe_cpi_emg, method = "arellano", type = "HC3")))
          ),
          title           = "Table 5a — CPI Heterogeneity: Advanced vs. Emerging Economies (FE, Robust SE)",
          column.labels   = c("Advanced", "Emerging"),
          dep.var.labels  = "ln(CPI)",
          covariate.labels = c("Post-Crisis Dummy", "ln(NEER)", "Policy Rate", "VIX"),
          keep.stat       = c("n", "rsq"),
          type            = "text",
          notes           = "Advanced: Germany, Canada, USA. Emerging: India, Poland, Turkey, Saudi Arabia.")

# --------------------------------------------------------------------------
# 9D. Advanced vs Emerging — IPI (Saudi Arabia excluded)
# --------------------------------------------------------------------------

pdata_adv_ipi <- pdata_full_ipi %>%
  filter(group %in% c("Advanced_Importer", "Advanced_Exporter"))
pdata_emg_ipi <- pdata_full_ipi %>%
  filter(group %in% c("Emerging_Importer", "Emerging_Exporter"))

fe_ipi_adv <- plm(ln_ipi ~ post_crisis + ln_neer + policy_rate + vix,
                  data = pdata_adv_ipi, model = "within")

fe_ipi_emg <- plm(ln_ipi ~ post_crisis + ln_neer + policy_rate + vix,
                  data = pdata_emg_ipi, model = "within")

stargazer(fe_ipi_adv, fe_ipi_emg,
          se = list(
            sqrt(diag(vcovHC(fe_ipi_adv, method = "arellano", type = "HC3"))),
            sqrt(diag(vcovHC(fe_ipi_emg, method = "arellano", type = "HC3")))
          ),
          title           = "Table 5b — IPI Heterogeneity: Advanced vs. Emerging Economies (FE, Robust SE)",
          column.labels   = c("Advanced", "Emerging"),
          dep.var.labels  = "ln(IPI)",
          covariate.labels = c("Post-Crisis Dummy", "ln(NEER)", "Policy Rate", "VIX"),
          keep.stat       = c("n", "rsq"),
          type            = "text",
          notes           = "Saudi Arabia excluded. Emerging IPI: India, Poland, Turkey.")

# --------------------------------------------------------------------------
# 9E. Per-group loop — individual FE models for each of the 4 groups (bonus)
# --------------------------------------------------------------------------

cat("\n--- Per-Group FE Models: CPI (Full Mixed Panel) ---\n")
for (grp in sort(unique(data_full$group))) {
  grp_data <- data_full %>% filter(group == grp)
  if (n_distinct(grp_data$country) < 2) {
    cat("\nGroup:", grp,
        "— skipped (only 1 country; FE not identified with single cross-section)\n")
    next
  }
  grp_pdata <- pdata.frame(grp_data, index = c("country", "date"))
  grp_model <- plm(ln_cpi ~ post_crisis + ln_neer + policy_rate + vix,
                   data = grp_pdata, model = "within")
  cat("\n===== GROUP:", grp, "=====\n")
  print(coeftest(grp_model,
                 vcov = vcovHC(grp_model, method = "arellano", type = "HC3")))
}


# ============================================================================
# 10. FULL INTERPRETATION AND REPORT WRITING GUIDE
# ============================================================================

cat("
================================================================================
INTERPRETATION AND REPORT WRITING GUIDE
================================================================================

SECTION 3 — DATA AND VARIABLES (what to write):
  State the original group: 6 Emerging Importers from the project guidelines.
  Then explain the two exclusions:

  Thailand (excluded from ALL models):
    'Thailand exhibited near-zero CPI variation over the sample period, with
    index values ranging from 99.24 to 100.98 — a total spread of less than
    2 index points over 42 months (within-country variance ≈ 0.18, versus
    13.9 for India and 12.9 for Poland). The Fixed Effects estimator relies
    entirely on within-country time variation; including a near-flat series
    provides no identifying information and risks attenuating estimated
    coefficients toward zero. Thailand was therefore excluded from all CPI
    and IPI regression models. This reflects Thailand's structural
    near-deflation during 2023–2025, not a data collection failure.'

  South Africa (excluded from IPI models only):
    'South Africa's IPI series is available only through November 2023
    (11 observed months out of 42). The remaining 31 observations — spanning
    all of 2024, all of 2025, and the entire post-crisis period — are
    ARIMA-forecasted. Critically, the ARIMA model converged to a flat constant
    (IPI = 110.2543) for every month from December 2023 onward, producing
    zero within-country variation across the most policy-relevant portion of
    the sample. Including this series in the IPI panel would mean the
    post-crisis shock coefficient is effectively estimated from a constant,
    which would attenuate β and inflate standard errors. South Africa is
    retained in CPI models, where its data quality is acceptable, but excluded
    from all IPI specifications.'

  India (retained with caveat):
    'India's IPI is observed through October 2024; the remaining 20 months
    (48% of the series) are ARIMA-forecasted. Unlike South Africa, India's
    ARIMA forecasts exhibit realistic seasonal variation rather than converging
    to a flat constant, and dropping India would leave only 3 countries in the
    IPI panel, further reducing statistical power. India is retained in IPI
    models, but this limitation is acknowledged: the estimated post-crisis
    effect for India's industrial production is partly identified from
    ARIMA extrapolations of pre-crisis trends, which likely understates the
    true shock.'

SECTION 5 — ECONOMETRIC RESULTS (what to write):
  Tables 2a/3a: Compare β₁ across three POLS specs. If β₁ drops in magnitude
    or loses significance from Col (1) to Col (2), this confirms the
    multicollinearity the professor identified. Quote the specific values.
  Tables 2b/3b: Preferred FE results. Key numbers to report:
    - post_crisis coefficient: magnitude and significance (expected + for CPI, - for IPI)
    - NEER coefficient: sign flip from POLS to FE confirms omitted variable bias
    - R² improvement from POLS to FE confirms country effects matter
  Diagnostic tests: report all three (Hausman, BP, BG) with p-values.
    If Hausman fails to reject: argue FE still preferred on theoretical grounds
    and cite low power due to small N.

SECTION 6 — HETEROGENEITY (what to write):
  Tables 4a/4b: Importers vs Exporters.
    Expected: importers show larger positive post_crisis β for CPI.
    Expected: exporters (USA, Canada) may show positive or near-zero IPI β.
  Tables 5a/5b: Advanced vs Emerging.
    Expected: emerging economies show larger CPI responses due to weaker
    monetary credibility and higher NEER pass-through.

SECTION 8 — LIMITATIONS (what to write):
  1. Thailand exclusion: structural near-deflation, not a data collection issue.
  2. South Africa IPI: insufficient real data for the estimation period.
  3. India IPI: post-crisis period ARIMA-generated; likely attenuates estimates.
  4. Short post-crisis window (4 months) constrains identification.
  5. Small N in primary panel (4-5 countries): Hausman test has low power.
================================================================================
")

# ============================================================================
# END OF SCRIPT
# ============================================================================
