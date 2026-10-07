# Data and code: Environment and host specificity shape the seasonal dynamics of co-occurring amphibian pathogens

Arik M. Hartmann, Max L. Maddox, Cory C. McKinstry, Steven M. Klioze, Steve A. Johnson and Ana V. Longo

Data and R code supporting the article "Environment and host specificity shape the seasonal dynamics of co-occurring amphibian pathogens", submitted to *Biology Letters*

- Archived copy (DOI): TBD
- Development copy: https://github.com/amphibiarik/Hartmann_et_al_2026_Biology_Letters
- Contact: Arik M. Hartmann (arikhartmann@gmail.com)

## Contents

| File | Format | Description |
|---|---|---|
| `Hartmann_BL_data_clean.csv` | Comma-separated values (UTF-8, one header row, 2,119 rows x 15 columns) | One row per sampled amphibian: survey date and wetland, 7-day weather summaries, host species and metamorphic stage, and mean qPCR estimates and status for *Batrachochytrium dendrobatidis* (Bd) and ranavirus (Rv). |
| `Hartmann_BioLetters_analysis_figures.R` | R script (plain text) | Reproduces every analysis, table and figure in the article and its electronic supplementary material. |
| `README.md` | Markdown | This description of the data and code. |

## How the data were collected

From December 2019 to March 2022 we sampled amphibian communities at five wetlands at Ordway-Swisher Biological Station (OSBS), north-central Florida, USA (29.6729° N, 82.0331° W). The wetlands are semi-permanent, fishless sinkhole ponds occupied by the striped newt (*Notophthalmus perstriatus*). Each wetland was surveyed approximately monthly when standing water was present (91 survey events on 85 dates).

Surveys combined dip-netting within 3 m of the pond margin (aquatic stages) with time-constrained searches of the wetland perimeter (terrestrial stages). Salamanders and post-metamorphic anurans were swabbed across the ventral surface, including the gular region and vent. Tadpoles were swabbed on the oral disc for Bd and on the vent for Rv. Swabs were stored individually in sterile 2 ml tubes at −20 °C. Animals were handled with clean gloves, and equipment was disinfected with 5% bleach between sites. Anurans were classed as post-metamorphic after complete tail resorption (Gosner stage 46); striped newts were classed from external morphology (gill presence, skin texture and secondary sexual characteristics).

DNA was extracted with PrepMan Ultra (Applied Biosystems). Bd and Rv were quantified in separate TaqMan qPCR assays (Boyle et al. 2004 for Bd; Allender et al. 2013 for Rv). Samples, standards and negative controls were run in duplicate, and samples that amplified in only one replicate were run a third time. Rv samples amplifying after 40 cycles were scored negative. Bd samples were scored positive when ITS copy number exceeded one. Loads are the mean copy number across replicate reactions, corrected for the sample dilution factor.

Weather data come from the OSBS automated weather station (https://osbs-met.ifas.ufl.edu/vdv.php/historical/46). For each survey we averaged daily values over the seven days preceding the survey.


## Column dictionary: `data/Hartmann_BL_data_clean.csv`

There are no missing values. Zeros in the two load columns mean the pathogen was not detected.

| Column | Type | Units / values | Description |
|---|---|---|---|
| `date` | date (text) | `M/D/YYYY`, e.g. `12/6/2019` | Survey date. |
| `pond_ID` | text | `AQ`, `BP`, `FP`, `OSP`, `PL` | Wetland code (see table above). A survey event is one wetland on one date. |
| `avg_temp_C` | numeric | °C | Mean of daily mean air temperature over the 7 days before the survey. |
| `precip_avg` | numeric | mm day⁻¹ | Mean daily precipitation over the 7 days before the survey. |
| `species` | text | 4-letter code (table below) | Host species. |
| `ID` | text | Number, or `OS` + number | Field identifier of the swabbed individual. Unique within a species except one striped newt (`NOPE`, `OS102`) sampled on two dates; the script builds a unique `individual_ID` from `spp` + `ID`. |
| `metam_stage` | text | `pre`, `post` | Metamorphic stage: `pre` = larva/tadpole (aquatic, pre-metamorphic); `post` = metamorphosed juvenile or adult. |
| `rv_ge` | numeric | viral genome copies (mean of replicate qPCR reactions) | Ranavirus load. 0 = not detected. Rounded to 2 decimal places. |
| `bd_its` | numeric | Bd ITS copies (mean of replicate qPCR reactions) | Bd load. 0 = not detected. Rounded to 2 decimal places. |
| `rv_status` | integer | 0 = negative, 1 = positive | Ranavirus infection status from qPCR. |
| `bd_status` | integer | 0 = negative, 1 = positive | Bd infection status from qPCR. |

### Species codes (`species`)

| Code | Species | n | Code | Species | n |
|---|---|---|---|---|---|
| ACGR | *Acris gryllus* | 470 | LICA | *Lithobates catesbeianus* | 112 |
| ANQU | *Anaxyrus quercicus* | 1 | LICP | *Lithobates capito* | 208 |
| ANTE | *Anaxyrus terrestris* | 9 | LIGR | *Lithobates grylio* | 27 |
| EUPL | *Eleutherodactylus planirostris* | 3 | LISP | *Lithobates sphenocephalus* | 167 |
| EUQU | *Eurycea quadridigitata* | 13 | NOPE | *Notophthalmus perstriatus* | 817 |
| GACA | *Gastrophryne carolinensis* | 3 | PSOC | *Pseudacris ocularis* | 7 |
| HYCI | *Dryophytes cinereus* | 5 | SCHO | *Scaphiopus holbrookii* | 3 |
| HYFE | *Dryophytes femoralis* | 110 | | | |
| HYGR | *Dryophytes gratiosus* | 162 | | | |
| HYSQ | *Dryophytes squirellus* | 2 | | | |

Abbreviations: Bd, *Batrachochytrium dendrobatidis*; Rv, ranavirus; ITS, internal transcribed spacer (Bd rDNA target); ge, genome equivalents; qPCR, quantitative PCR; OSBS, Ordway-Swisher Biological Station.

## How the data were processed and analysed

All analysis and visualization steps are in `R/Hartmann_BioLetters_analysis_figures.R`, in the order of the numbered sections in the script.

1. **Data preparation.** Dates are parsed; day of year, year and month are derived; a coinfection indicator (Bd and Rv both positive) and a unique individual ID are created. Analyses of infection probability use the seven species with ≥100 sampled individuals (ACGR, HYFE, HYGR, LICA, LICP, LISP, NOPE; 2,046 individuals). Temperature (`avg_temp_C`) and precipitation (`precip_avg`) are z-standardised. Loads are log10(x + 1) transformed.
2. **Summary statistics (Table S1).** Prevalence of Bd, Rv and coinfection by species for all 17 species, with exact binomial 95% confidence intervals.
3. **Infection probability.** Separate binomial GLMMs (glmmTMB) for Bd and Rv status with species, metamorphic stage, 7-day mean temperature and 7-day mean precipitation as fixed effects and survey event as a random intercept. Type II Wald χ² tests (car), Tukey-adjusted pairwise comparisons of estimated marginal means (emmeans), and residual diagnostics (DHARMa: uniformity, dispersion, outliers; performance: singularity, variance inflation factors). Species × stage and temperature × precipitation interactions were tested and not retained.
4. **Pathogen load.** Gaussian GLMMs of log10 load among qPCR-positive individuals (species with ≥5 positives), with the same fixed and random effects and the same tests and diagnostics.
5. **Seasonal timing and strength.** A joint binomial GLMM of Bd and Rv status with pathogen-specific annual harmonics (sine and cosine of day of year) and pathogen-specific effects of species, stage and sampling season (July–June), with random intercepts for wetland and for survey event within pathogen. Harmonic order was chosen by AIC. A likelihood-ratio test compares this model with one sharing a single seasonal cycle. Peak dates, the Rv − Bd lag and seasonal amplitudes, with 95% intervals, come from 10,000 draws of the fixed effects (MASS::mvrnorm). *D. femoralis* is excluded because it was never Bd-positive. Also: within-species lags for ACGR and NOPE, prevalence by sampling season, and a Spearman correlation of survey-level Bd and Rv prevalence (surveys with ≥10 animals).
6. **Coinfection.** For *A. gryllus* and *N. perstriatus*, a stratified permutation test: Rv status is permuted among individuals of the same species, survey event and life stage, with Bd status held fixed (10,000 permutations; one-sided P = proportion of permutations with at least as many coinfections as observed).
7. **Figures.** Figure 1 (A–D), Figure 2 (A–D) and supplementary Figures S1–S2 (ggplot2, patchwork, ggdist). `ggsave()` lines are commented out; uncomment them to write TIFF files.

Random seeds are set (`set.seed(123)`) before every simulation and permutation step, so results reproduce exactly.

## Software

Analyses were run in R 4.4.3 (R Core Team 2025) on [operating system].

| Package | Used for |
|---|---|
| tidyverse (dplyr, tidyr, ggplot2, purrr, tibble) | Data handling and plotting |
| lubridate | Dates |
| glmmTMB | [x.y.z] | Mixed models |
| car | Type II Wald tests |
| emmeans | Marginal means and contrasts |
| DHARMa | Residual diagnostics |
| performance | Singularity and collinearity checks |
| broom.mixed | Model tables |
| patchwork | Figure layout |
| MASS | Multivariate normal draws |
| ggdist | Half-eye plots (Fig. 2D) |
| scales | Axis labels |

To reproduce: open R, install the packages above, and run Hartmann_BioLetters_analysis_figures.R. The permutation and simulation steps take a few minutes.

## Licence and access

## References

- Allender MC, Bunick D, Mitchell MA. 2013 Development and validation of TaqMan quantitative PCR for detection of frog virus 3-like virus in eastern box turtles (*Terrapene carolina carolina*). *J. Virol. Methods* 188, 121–125.
- Boyle DG, Boyle DB, Olsen V, Morgan JAT, Hyatt AD. 2004 Rapid quantitative detection of chytridiomycosis (*Batrachochytrium dendrobatidis*) in amphibian samples using real-time Taqman PCR assay. *Dis. Aquat. Org.* 60, 141–148.
- R Core Team. 2025 R: a language and environment for statistical computing. Vienna, Austria: R Foundation for Statistical Computing.
