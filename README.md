
# The Geopolitics of Oil: Macroeconomic Impact of the 2026 Gulf Crisis

Advanced Econometrics project (panel data analysis in R)

**Team:** Maycem Azaza · Rahma Haddouchi · Noursine Amira · Lina Brahmi
**Period covered:** January 2023 – June 2026 (42 monthly observations per country)

---

## Overview

On 27 February 2026, escalating tensions in the Persian Gulf triggered a sharp oil price shock. This project treats the event as a quasi-natural experiment and estimates its effect on:

- **Inflation** (log CPI)
- **Industrial production** (log IPI)

The primary sample is a panel of emerging oil-importing economies. A broader mixed panel is used to study heterogeneity (importers vs. exporters, advanced vs. emerging).

## Main Findings

| Outcome | Preferred model (FE) | Effect of the crisis |
|---|---|---|
| CPI (5 emerging importers) | Post-crisis dummy +0.031 (significant) | about +3.1% |
| IPI (4 emerging importers) | Post-crisis dummy −0.022 (significant) | about −2.2% |

- The pattern is consistent with a supply-side **stagflationary shock**.
- The nominal effective exchange rate (NEER) acts as a shock absorber: currency appreciation lowers CPI and supports output.
- Emerging importers are hit harder than advanced economies and oil exporters.
- Fixed Effects is preferred over POLS and RE (Breusch-Pagan rejects POLS; Hausman rejects RE for CPI).

## Repository Contents

| File | Description |
|---|---|
| `Rcode.R` | Full analysis script: data prep, descriptives, figures, models, diagnostics, heterogeneity |
| `Group_Dataset.csv` | Emerging Importers panel (6 countries, 252 rows); used for all primary models |
| `Full mixed dataset.csv` | Seven-country mixed panel (294 rows); used only for heterogeneity analysis |
| Project report (PDF) | Full write-up of methodology, results and policy implications |

> **Note on file names:** the script reads `panel_dataset_FINAL_merged_fixed__2_.csv`. The file provided here is called `Full mixed dataset.csv`. Either rename the CSV to match, or change the `read.csv()` call in section 2B of the script.

## Data

Monthly data. Sources listed in the report: IMF IFS, OECD, FRED, EIA, BIS, national central banks.

| Variable | Description |
|---|---|
| `cpi` | Consumer Price Index (log-transformed) |
| `ipi` | Industrial Production Index (log-transformed) |
| `brent` | Brent crude price; `oil_pct_change` is its monthly % change |
| `neer` | Nominal Effective Exchange Rate (log-transformed) |
| `policy_rate` | Central bank policy rate (%) |
| `vix` | CBOE Volatility Index (global risk proxy) |
| `post_crisis` | Dummy = 1 from March 2026 onward |
| `group`, `oil_balance` | Country classification (advanced/emerging, importer/exporter) |
| `*_flag` columns | Mark observed vs. ARIMA-forecasted values |
| `ipi_model_eligible` | `included` / `excluded` flag for the heterogeneity IPI panel (Saudi Arabia excluded) |

Values that were not yet available at collection time (notably April–June 2026) were filled by country-level ARIMA forecasts.

## Sample Decisions

| Country | CPI models | IPI models | Reason |
|---|---|---|---|
| India | Included | Included* | *48% of IPI is ARIMA-forecasted; limitation acknowledged |
| Indonesia | Included | Included | Data quality acceptable |
| Poland | Included | Included | Data quality acceptable |
| South Africa | Included | **Excluded** | 74% of IPI is ARIMA, flat at 110.25 from Dec 2023 |
| Thailand | **Excluded** | **Excluded** | Near-zero CPI variance (range 99.24–100.98) |
| Turkey | Included | Included | Data quality acceptable |

Final panels: **CPI** = 5 countries (210 obs); **IPI** = 4 countries (168 obs).

## Methodology

### Identification strategy

Brent prices are identical across countries in any month, so `oil_pct_change` is collinear with time effects and cannot be identified under the FE/RE within transformation. Two tracks are used:

- **Pooled OLS:** includes `oil_pct_change` directly (direct price pass-through).
- **FE and RE:** use the `post_crisis` dummy as a reduced-form proxy for the shock.
- **Robustness:** a POLS model with both variables shows the multicollinearity (the dummy absorbs the oil coefficient).

### Models (for each of CPI and IPI)

1. POLS, oil only
2. POLS, oil + crisis dummy (robustness)
3. POLS, crisis dummy only
4. Fixed Effects (preferred)
5. Random Effects

Controls: `ln(NEER)`, `policy_rate`, `VIX`.

### Diagnostics and standard errors

- Hausman test (FE vs. RE)
- Breusch-Pagan LM test (POLS vs. RE)
- Breusch-Godfrey test (serial correlation)
- Breusch-Pagan test (heteroskedasticity)
- Robust SEs: HC3 for POLS, Arellano HC3 for FE/RE

### Heterogeneity analysis

FE models estimated separately on the mixed panel for:

- Oil importers vs. exporters
- Advanced vs. emerging economies
- Individual groups (bonus loop)

## How to Run

### Requirements

R with the following packages:

```r
install.packages(c("tidyverse", "plm", "lmtest", "stargazer",
                   "ggplot2", "car", "sandwich"))
```

### Steps

1. Place `Rcode.R`, `Group_Dataset.csv` and the mixed-panel CSV in the same folder.
2. Make sure the mixed-panel file name matches the one in the script (see note above).
3. Set the working directory to that folder, e.g. `setwd("path/to/project")`.
4. Run the script:

```r
source("Rcode.R")
```

The script prints all tables to the console (`stargazer`, text format), runs the diagnostics, and saves the figures.

### Outputs

Nine figures are saved as PNG:

| Figure | Content |
|---|---|
| `fig01_brent_price.png` | Brent crude price with crisis onset |
| `fig02_lncpi_all6.png` | Log CPI, all 6 emerging importers |
| `fig03_lncpi_sample.png` | Log CPI, CPI estimation sample |
| `fig04_lnipi_sample.png` | Log IPI, IPI estimation sample |
| `fig05_sa_ipi_arima.png` | South Africa IPI: observed vs. ARIMA |
| `fig06_lncpi_bygroup.png` | Average log CPI by group |
| `fig07_lnipi_bygroup.png` | Average log IPI by group |
| `fig08_cpi_brent_divergence.png` | CPI vs. Brent (dual axis) |
| `fig09_oil_shock_monthly.png` | Monthly % change in Brent |

## Limitations

- Thailand is excluded because of near-flat CPI (a structural feature, not a data error).
- South Africa's IPI is largely unusable (flat ARIMA values).
- India's post-crisis IPI is ARIMA-extrapolated and likely biases its effect toward zero.
- The post-crisis window is short (March–June 2026), so second-round effects are probably under-captured.
- Small cross-section (N = 4–5) gives the Hausman test low power and limits generalisation.
