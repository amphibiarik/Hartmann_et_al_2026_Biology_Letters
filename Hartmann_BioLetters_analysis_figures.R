# ============================================================
# Hartmann et al. submitted to Biology Letters
# Bd and Ranavirus co-emergence dynamics in amphibian communities

#  1. Data preparation
#  2. Summary statistics (Table S1)
#  3. Prevalence models (species, stage, temperature, precipitation)
#  4. Load models (qPCR-positive individuals)
#  5. Seasonal timing and strength (joint Bd-Rv harmonic GLMM)
#  6. Coinfection permutation test
#  7. Figures: Fig. 1 (A-D), Fig. 2 (A-D), Fig. S1, Fig. S2
# ============================================================

# Packages ----
library(tidyverse)
library(lubridate)
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(broom.mixed)
library(performance)
library(patchwork)

# Functions called with :: so they do not mask dplyr: MASS (mvrnorm), ggdist (stat_halfeye)


# 1. Data preparation ----
# Variables used: date, pond_ID, avg_temp_C, precip_avg, species, ID, metam_stage,
# bd_its (Bd load) mean from duplicate qPCR runs, rv_ge (Rv load) mean from duplicate qPCR runs
# bd_status, rv_status, = 1 if qPCR mean was above detectable thresholds 


d_raw <- read.csv("Hartmann_BL_data_clean.csv") %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"),
         DOY = yday(date), year = year(date), month = month(date),
         site = factor(pond_ID), species = factor(species), stage = factor(metam_stage),
         coinfected = case_when(bd_status == 1 & rv_status == 1 ~ 1,
                                !is.na(bd_status) & !is.na(rv_status) ~ 0)) %>%
  # Unique IDs = species code + swab series (duplicates get a suffix)
  group_by(species, ID) %>%
  mutate(individual_ID = if (n() > 1) paste(species, ID, row_number(), sep = "_")
         else paste(species, ID, sep = "_")) %>%
  ungroup()

d_raw
# 2,119 individuals, 17 species

# Species with adequate sample sizes (Eurycea: only 2 wetlands and 6 surveys, excluded)
species_keep <- c("ACGR", "HYFE", "HYGR", "LICA", "LICP", "LISP", "NOPE")

d <- d_raw %>%
  mutate(survey_event = interaction(site, date, drop = TRUE),
         log_bd_load = log10(bd_its + 1),
         log_rv_load = log10(rv_ge + 1)) %>%
  filter(species %in% species_keep) %>%
  mutate(species = droplevels(species), stage = droplevels(stage),
         avg_temp_z = as.numeric(scale(avg_temp_C)),
         precip_z = as.numeric(scale(precip_avg)))

#Filtered data 2,046 individuals, 7 species

prev_base   <- d %>% filter(complete.cases(species, stage, avg_temp_z, precip_z, survey_event))
bd_prev_dat <- prev_base %>% filter(!is.na(bd_status)) %>% droplevels()
rv_prev_dat <- prev_base %>% filter(!is.na(rv_status)) %>% droplevels()

d %>% summarise(n = n(), n_individuals = n_distinct(individual_ID), n_species = n_distinct(species),
                bd_positive = sum(bd_status == 1, na.rm = TRUE), rv_positive = sum(rv_status == 1, na.rm = TRUE))

# 2. Summary statistics ----
binom_ci <- function(x, n) binom.test(x, n)$conf.int

overall_prev <- tibble(
  pathogen = c("Bd", "Rv", "Co-infected"),
  positive = c(sum(d_raw$bd_status == 1, na.rm = TRUE), sum(d_raw$rv_status == 1, na.rm = TRUE),
               sum(d_raw$coinfected == 1, na.rm = TRUE)),
  tested   = c(sum(!is.na(d_raw$bd_status)), sum(!is.na(d_raw$rv_status)), sum(!is.na(d_raw$coinfected)))) %>%
  rowwise() %>%
  mutate(prevalence = positive / tested,
         lower = binom_ci(positive, tested)[1], upper = binom_ci(positive, tested)[2]) %>%
  ungroup()
overall_prev

# Table S1: prevalence by species (all 17 species)
fmt_prev <- function(x, n) { ci <- binom_ci(x, n); sprintf("%.1f%% (%.1f–%.1f)", 100 * x / n, 100 * ci[1], 100 * ci[2]) }

species_table <- d_raw %>%
  group_by(species) %>%
  summarise(N_sampled = n(),
            Bd_pos = sum(bd_status == 1, na.rm = TRUE), Bd_n = sum(!is.na(bd_status)),
            Rv_pos = sum(rv_status == 1, na.rm = TRUE), Rv_n = sum(!is.na(rv_status)),
            Co_pos = sum(coinfected == 1, na.rm = TRUE), Co_n = sum(!is.na(coinfected)), .groups = "drop") %>%
  rowwise() %>%
  transmute(Species = species, `N sampled` = N_sampled,
            `N Bd positive` = Bd_pos, `Bd prevalence (95% CI)` = fmt_prev(Bd_pos, Bd_n),
            `N Rv positive` = Rv_pos, `Rv prevalence (95% CI)` = fmt_prev(Rv_pos, Rv_n),
            `N co-infected` = Co_pos, `Co-infection prevalence (95% CI)` = fmt_prev(Co_pos, Co_n)) %>%
  ungroup()

#Inspect, uncomment below to save summary table
species_table
#write.csv(species_table, "species_pathogen_summary.csv", row.names = FALSE)

# Model diagnostics helper for DHARMa residuals + performance checks on glms
check_glmm <- function(model, seed = 123, nsim = 2000) {
  set.seed(seed)
  res <- DHARMa::simulateResiduals(model, n = nsim)
  plot(res)
  list(singularity = performance::check_singularity(model),
       collinearity = performance::check_collinearity(model),
       uniformity = DHARMa::testUniformity(res), dispersion = DHARMa::testDispersion(res),
       outliers = DHARMa::testOutliers(res), random_effects = VarCorr(model))
}

# 3. Prevalence models ----
bd_prev_final <- glmmTMB(bd_status ~ species + stage + avg_temp_z + precip_z + (1 | survey_event),
                         family = binomial, data = bd_prev_dat)
rv_prev_final <- glmmTMB(rv_status ~ species + stage + avg_temp_z + precip_z + (1 | survey_event),
                         family = binomial, data = rv_prev_dat)

car::Anova(bd_prev_final, type = 2)  # species, stage, temperature, precipitation
car::Anova(rv_prev_final, type = 2)  # species, temperature

broom.mixed::tidy(bd_prev_final, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)  # HYFE large: never Bd+
broom.mixed::tidy(rv_prev_final, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)

check_glmm(bd_prev_final)
check_glmm(rv_prev_final)  # outlier test significant, expected given Rv outbreak data

# Species and stage post hocs
emm_bd_species <- emmeans(bd_prev_final, ~ species, type = "response", weights = "proportional")
emm_rv_species <- emmeans(rv_prev_final, ~ species, type = "response", weights = "proportional")
pairs(emm_bd_species, adjust = "tukey")  # ACGR and NOPE higher; HYFE unstable
pairs(emm_rv_species, adjust = "tukey")  # LICP higher

emm_bd_stage <- emmeans(bd_prev_final, ~ stage, type = "response")
emm_rv_stage <- emmeans(rv_prev_final, ~ stage, type = "response")
pairs(emm_bd_stage)  # post-metamorphic OR = 1.73
pairs(emm_rv_stage)  # NS

# Optional species x stage and temperature x precipitation interactions (none significant)
# Uncomment below to re-examine 
#bd_prev_int <- update(bd_prev_final, . ~ species * stage + avg_temp_z * precip_z + (1 | survey_event))
#rv_prev_int <- update(rv_prev_final, . ~ species * stage + avg_temp_z * precip_z + (1 | survey_event))
#car::Anova(bd_prev_int, type = 2)
#car::Anova(rv_prev_int, type = 2)
#anova(bd_prev_final, bd_prev_int)
#anova(rv_prev_final, rv_prev_int)  # better fit, but no interaction term significant

# 4. Load models (positive individuals; species with >= 5 positives) ----
bdL <- d %>% filter(bd_status == 1, bd_its > 0, !is.na(bd_its)) %>%
  group_by(species) %>% filter(n() >= 5) %>% ungroup() %>% mutate(species = droplevels(species))

rvL <- d %>% filter(rv_status == 1, rv_ge > 0, !is.na(rv_ge)) %>%
  group_by(species) %>% filter(n() >= 5) %>% ungroup() %>% mutate(species = droplevels(species))

bdL %>% count(species)
rvL %>% count(species)

bdL_final <- glmmTMB(log_bd_load ~ species + stage + avg_temp_z + precip_z + (1 | survey_event),
                     family = gaussian, data = bdL)

rvL_final <- glmmTMB(log_rv_load ~ species + stage + avg_temp_z + precip_z + (1 | survey_event),
                     family = gaussian, data = rvL)

car::Anova(bdL_final, type = 2)
car::Anova(rvL_final, type = 2)
broom.mixed::tidy(bdL_final, effects = "fixed", conf.int = TRUE)
broom.mixed::tidy(rvL_final, effects = "fixed", conf.int = TRUE)

check_glmm(bdL_final)
check_glmm(rvL_final)  # outlier/deviation tests significant, expected with Rv outbreak loads during disease events

pairs(emmeans(bdL_final, ~ species), adjust = "tukey")  # none significant
pairs(emmeans(rvL_final, ~ species), adjust = "tukey")  # LICP and NOPE
pairs(emmeans(bdL_final, ~ stage))  # post-metamorphic higher
pairs(emmeans(rvL_final, ~ stage))  # NS

# Post hoc tables (raw CSVs + formatted Table S2) ----
# All species and stage pairs() from sections 3-4. Prevalence contrasts are odds ratios;
# load contrasts are differences in log10 load. Species CIs are Tukey-adjusted, like the p-values.
# HYFE was never Bd+, so its Bd prevalence contrasts are not estimable (complete separation):
# they stay in the raw CSVs with separation = TRUE and are shown as "n.e." in Table S2.
posthoc <- list(
  bd_prev_species = pairs(emm_bd_species, adjust = "tukey"),
  rv_prev_species = pairs(emm_rv_species, adjust = "tukey"),
  bd_prev_stage   = pairs(emm_bd_stage),
  rv_prev_stage   = pairs(emm_rv_stage),
  bd_load_species = pairs(emmeans(bdL_final, ~ species), adjust = "tukey"),
  rv_load_species = pairs(emmeans(rvL_final, ~ species), adjust = "tukey"),
  bd_load_stage   = pairs(emmeans(bdL_final, ~ stage)),
  rv_load_stage   = pairs(emmeans(rvL_final, ~ stage)))

# One tidy tibble per table, with the same column names whatever the model family
tidy_pairs <- function(x) {
  summary(x, infer = TRUE) %>% as.data.frame() %>% as_tibble() %>%
    rename_with(~ case_when(.x %in% c("odds.ratio", "ratio", "estimate") ~ "estimate",
                            .x %in% c("asymp.LCL", "lower.CL") ~ "lower",
                            .x %in% c("asymp.UCL", "upper.CL") ~ "upper",
                            .x %in% c("z.ratio", "t.ratio") ~ "statistic",
                            TRUE ~ .x)) %>%
    mutate(contrast = as.character(contrast))
}

posthoc_raw <- map(posthoc, tidy_pairs) %>%
  imap(~ mutate(.x, separation = str_starts(.y, "bd_prev") & str_detect(contrast, "HYFE")))

posthoc_all <- bind_rows(posthoc_raw, .id = "table") %>%
  separate_wider_delim(table, "_", names = c("pathogen", "response", "term")) %>%
  mutate(pathogen = str_to_title(pathogen),
         response = if_else(response == "prev", "Infection probability", "log10 load"),
         term = str_to_title(term))
posthoc_all %>% print(n = Inf)

#Uncomment below to save one raw CSV per post hoc table (plus all tables stacked)
#iwalk(posthoc_raw, ~ write.csv(.x, paste0("posthoc_", .y, ".csv"), row.names = FALSE))
#write.csv(posthoc_all, "posthoc_all_raw.csv", row.names = FALSE)

# Tables S2-S9: one formatted table per post hoc, for the supplement (species names as in Table S1)
species_names <- c(ACGR = "Acris gryllus", HYFE = "Dryophytes femoralis", HYGR = "Dryophytes gratiosus",
                   LICA = "Lithobates catesbeianus", LICP = "Lithobates capito",
                   LISP = "Lithobates sphenocephalus", NOPE = "Notophthalmus perstriatus")
fmt_p <- function(p) if_else(p < 0.001, "<0.001", sprintf("%.3f", p))

posthoc_info <- tibble(
  table = names(posthoc),
  number = paste0("S", seq_along(posthoc) + 1),
  caption = c("Pairwise comparisons of Bd infection probability among species.",
              "Pairwise comparisons of Rv infection probability among species.",
              "Comparison of Bd infection probability between post- and pre-metamorphic individuals.",
              "Comparison of Rv infection probability between post- and pre-metamorphic individuals.",
              "Pairwise comparisons of Bd load among species (Bd-positive individuals).",
              "Pairwise comparisons of Rv load among species (Rv-positive individuals).",
              "Comparison of Bd load between post- and pre-metamorphic individuals (Bd-positive individuals).",
              "Comparison of Rv load between post- and pre-metamorphic individuals (Rv-positive individuals)."),
  footnote = case_when(
    str_detect(table, "prev") ~ "Odds ratios from the binomial GLMM, on the response scale; z tests.",
    TRUE ~ "Differences in log10(load + 1) from the Gaussian GLMM; t tests.") %>%
    paste(if_else(str_detect(table, "species"), "P-values and 95% CIs are Tukey-adjusted.", "")) %>%
    paste(if_else(table == "bd_prev_species",
                  "n.e. = not estimable: D. femoralis was never Bd positive (complete separation).", "")) %>%
    str_squish())

fmt_posthoc <- function(x, table) {
  prev <- str_detect(table, "prev")
  x %>%
    mutate(Contrast = str_replace_all(contrast, species_names) %>% str_replace_all(" [/-] ", " vs. "),
           est = if_else(separation, "n.e.",
                         sprintf(if (prev) "%.2f (%.2f–%.2f)" else "%.2f (%.2f to %.2f)", estimate, lower, upper)),
           stat = if_else(separation, "–", sprintf("%.2f", statistic)),
           p = if_else(separation, "–", fmt_p(p.value))) %>%
    select(Contrast, est, stat, p) %>%
    set_names(c("Contrast", if (prev) c("Odds ratio (95% CI)", "z") else c("Difference (95% CI)", "t"), "p"))
}

posthoc_fmt <- imap(posthoc_raw, fmt_posthoc)
posthoc_fmt

#Uncomment below to save Tables S2-S9 as one CSV each and as one Word file (.docx; needs flextable and officer)
pwalk(posthoc_info, function(table, number, ...)
  write.csv(posthoc_fmt[[table]], paste0("Table_", number, "_", table, ".csv"), row.names = FALSE, fileEncoding = "UTF-8"))
posthoc_ft <- function(x, footnote) {
  flextable::flextable(x) %>%
    flextable::bold(part = "header") %>%
    flextable::add_footer_lines(footnote) %>%
    flextable::font(fontname = "Times New Roman", part = "all") %>%
    flextable::fontsize(size = 10, part = "all") %>%
    flextable::padding(padding.top = 1, padding.bottom = 1, part = "all") %>%
    flextable::width(width = c(3.4, 1.7, 0.6, 0.7))  # inches; fits a 6.5 in text width
}
posthoc_docx <- officer::read_docx()
for (i in seq_len(nrow(posthoc_info))) {
  posthoc_docx <- posthoc_docx %>%
    officer::body_add_par(paste0("Table ", posthoc_info$number[i], ". ", posthoc_info$caption[i])) %>%
    flextable::body_add_flextable(posthoc_ft(posthoc_fmt[[posthoc_info$table[i]]], posthoc_info$footnote[i])) %>%
    officer::body_add_par("")
}
print(posthoc_docx, target = "Tables_S2-S9_posthoc.docx")

# 5. Seasonal timing and strength ----
# Joint binomial GLMM of Bd and Rv infection with pathogen-specific annual harmonics
# (sin/cos of day-of-year); LRT vs. one shared cycle; peak dates, lag and amplitude
# with 95% CIs from 10,000 draws of the fixed effects. Temperature and precipitation are
# omitted on purpose: this model describes WHEN infection peaks, not its correlates.
# HYFE excluded: never Bd-positive (complete separation in pathogen x species terms).
# season_yr = sampling season (Jul-Jun), so each winter stays in one level.
P <- 365.25

async_dat <- prev_base %>%
  filter(!is.na(bd_status), !is.na(rv_status), species != "HYFE") %>%
  droplevels() %>%
  mutate(season_yr = factor(ifelse(month >= 7, year, year - 1)),
         s1 = sin(2 * pi * DOY / P), c1 = cos(2 * pi * DOY / P),
         s2 = sin(4 * pi * DOY / P), c2 = cos(4 * pi * DOY / P)) %>%
  pivot_longer(c(bd_status, rv_status), names_to = "pathogen", values_to = "infected") %>%
  mutate(pathogen = factor(ifelse(pathogen == "bd_status", "Bd", "Rv"), levels = c("Bd", "Rv")),
         survey_path = interaction(survey_event, pathogen, drop = TRUE))

# Random effects: wetland, and survey event separately for each pathogen.
# (1 | individual_ID) had ~0 variance (AIC +2) so we omitted.
fit_async <- function(dat, seas, host = "species + stage + season_yr", shared = FALSE) {
  fixed <- if (shared) paste0("pathogen * (", host, ") + ", seas)
  else paste0("pathogen * (", host, " + ", seas, ")")
  glmmTMB(as.formula(paste0("infected ~ ", fixed, " + (1 | site) + (1 | survey_path)")),
          family = binomial, data = dat)
}

# Harmonic order: 1 (one peak per year) vs 2 (asymmetric/bimodal)
m_h1 <- fit_async(async_dat, "s1 + c1")
m_h2 <- fit_async(async_dat, "s1 + c1 + s2 + c2")
AIC(m_h1, m_h2)  # order 1: 2523.8 vs 2527.5
use_h2     <- isTRUE(AIC(m_h2) < AIC(m_h1))
m_full     <- if (use_h2) m_h2 else m_h1
seas_terms <- if (use_h2) c("s1", "c1", "s2", "c2") else c("s1", "c1")

# Pathogen-specific vs shared seasonal cycle
m_shared <- fit_async(async_dat, paste(seas_terms, collapse = " + "), shared = TRUE)
anova(m_shared, m_full)  # chi2 = 17.17, df = 2, P < 0.001

# Seasonal curves (logit scale) for a coefficient vector
doy_grid <- 1:365
H <- cbind(s1 = sin(2 * pi * doy_grid / P), c1 = cos(2 * pi * doy_grid / P),
           s2 = sin(4 * pi * doy_grid / P), c2 = cos(4 * pi * doy_grid / P))[, seas_terms, drop = FALSE]
wrap <- function(x) ((x + P / 2) %% P) - P / 2
seasonal_curves <- function(beta) {
  bd <- beta[seas_terms]
  rv <- bd + beta[paste0("pathogenRv:", seas_terms)]
  list(Bd = drop(H %*% bd), Rv = drop(H %*% rv))
}
peak_of <- function(curve) doy_grid[which.max(curve)]

summarise_peaks <- function(model, nsim = 10000, seed = 123) {
  b <- fixef(model)$cond
  set.seed(seed)
  sims <- MASS::mvrnorm(nsim, b, vcov(model)$cond)
  hat <- seasonal_curves(b)
  pk_hat <- c(Bd = peak_of(hat$Bd), Rv = peak_of(hat$Rv))
  pk <- t(apply(sims, 1, function(x) { cc <- seasonal_curves(x); c(peak_of(cc$Bd), peak_of(cc$Rv)) }))
  lag <- wrap(pk[, 2] - pk[, 1])
  # amplitude = half the peak-to-trough range of the seasonal component (logit scale)
  amp_of <- function(x) { cc <- seasonal_curves(x); c(diff(range(cc$Bd)), diff(range(cc$Rv))) / 2 }
  am_hat <- amp_of(b); am <- t(apply(sims, 1, amp_of)); am_ratio <- am[, 1] / am[, 2]
  ci <- function(v, h) (quantile(wrap(v - h), c(.025, .975)) + h) %% P  # circular-safe
  to_date <- function(x) format(as.Date(x - 1, origin = "2021-01-01"), "%d %b")
  tibble(
    quantity = c("Bd peak (DOY)", "Rv peak (DOY)", "Rv - Bd lag (days)",
                 "Bd amplitude (logit)", "Rv amplitude (logit)", "Bd / Rv amplitude ratio"),
    estimate = c(pk_hat["Bd"], pk_hat["Rv"], wrap(pk_hat["Rv"] - pk_hat["Bd"]), am_hat, am_hat[1] / am_hat[2]),
    lower = c(ci(pk[, 1], pk_hat["Bd"])[1], ci(pk[, 2], pk_hat["Rv"])[1], quantile(lag, .025),
              quantile(am[, 1], .025), quantile(am[, 2], .025), quantile(am_ratio, .025)),
    upper = c(ci(pk[, 1], pk_hat["Bd"])[2], ci(pk[, 2], pk_hat["Rv"])[2], quantile(lag, .975),
              quantile(am[, 1], .975), quantile(am[, 2], .975), quantile(am_ratio, .975)),
    p_rv_later = c(NA, NA, mean(lag > 0), NA, NA, NA)) %>%
    mutate(date_est = ifelse(grepl("peak", quantity), to_date(estimate), NA),
           date_lwr = ifelse(grepl("peak", quantity), to_date(lower), NA),
           date_upr = ifelse(grepl("peak", quantity), to_date(upper), NA))
}

async_peaks <- summarise_peaks(m_full)
async_peaks

# Bd peak 24 Feb (7 Feb-13 Mar); Rv peak 8 Feb (31 Dec - 19 Mar); lag -16 d (-58 to +27): timing NOT different, kinda expected
# Amplitude: Bd 2.73 (2.02-3.50), Rv 0.987 (0.56-1.51); ratio 2.76 (1.63-5.04): cycles differ in STRENGTH!! Very cool

# Within-species lags, separating pathogen timing from seasonal host turnover
# For ACGR and NOPE because they are most infected
within_species <- function(sp) {
  dat <- async_dat %>% filter(species == sp) %>% droplevels()
  m <- fit_async(dat, paste(seas_terms, collapse = " + "), host = "stage + season_yr")
  summarise_peaks(m) %>% mutate(species = sp, .before = 1)
}

bind_rows(within_species("ACGR"), within_species("NOPE"))

# ACGR lag +6 d (-43 to +54); NOPE lag -69 d (-128 to +13)

# DHARMa diagnostics
set.seed(123)
async_res <- DHARMa::simulateResiduals(m_full, n = 1000)
plot(async_res)
DHARMa::testUniformity(async_res); DHARMa::testDispersion(async_res)
m_full$sdr$pdHess  # TRUE
VarCorr(m_full)

# Interannual prevalence by sampling season (2019 = Dec 2019-Jun 2020 only)
prev_base %>%
  mutate(season_yr = ifelse(month >= 7, year, year - 1)) %>%
  group_by(season_yr) %>%
  summarise(n = n(), Bd = mean(bd_status), Rv = mean(rv_status), .groups = "drop")

# Bd 7.5%, 15.7%, 19.1%; Rv 16.7%, 27.4%, 9.99%

# Survey-level covariation of Bd and Rv prevalence (surveys with >= 10 animals)
survey_prev <- prev_base %>%
  group_by(survey_event) %>%
  summarise(n = n(), Bd = mean(bd_status), Rv = mean(rv_status), .groups = "drop") %>%
  filter(n >= 10)
cor.test(survey_prev$Bd, survey_prev$Rv, method = "spearman", exact = FALSE)  # rho = 0.28, P = 0.011

# 6. Coinfection permutation test ----
d %>% group_by(species) %>%
  summarise(n = n(), observed = sum(bd_status == 1 & rv_status == 1),
            expected = n * mean(bd_status) * mean(rv_status))

# 62 of 65 coinfections are in ACGR and NOPE!

coinf_dat <- d %>% filter(!is.na(bd_status), !is.na(rv_status), !is.na(stage))

# Permute Rv status within survey event x stage (Bd fixed), so each stratum keeps its Bd and Rv prevalence
coinf_perm <- function(species, nsim = 10000) {
  x <- coinf_dat[coinf_dat$species == species, ]
  obs <- sum(x$bd_status == 1 & x$rv_status == 1)
  sims <- replicate(nsim, {
    xp <- x %>% group_by(survey_event, stage) %>% mutate(rv_perm = sample(rv_status)) %>% ungroup()
    sum(xp$bd_status == 1 & xp$rv_perm == 1)
  })
  tibble(species = species, n = nrow(x), observed = obs, expected = mean(sims),
         lower = unname(quantile(sims, 0.025)), upper = unname(quantile(sims, 0.975)),
         p_upper = mean(sims >= obs))
}

set.seed(123)

coinf_perm_results <- bind_rows(coinf_perm("ACGR"), coinf_perm("NOPE"))
coinf_perm_results  # ACGR 37 vs 31.6, P = 0.028; NOPE 25 vs 24.8, P = 0.549

# 7. Figures ----
path_pal    <- c("Bd" = "#2F8AC4", "Rv" = "#DAA51B")
path_labels <- c("Bd" = expression(italic("B. dendrobatidis")), "Rv" = "Ranavirus")
spp_names <- c("ACGR" = "A. gryllus", "HYFE" = "D. femoralis", "HYGR" = "D. gratiosus",
               "EUQU" = "E. quadridigitata", "LICP" = "L. capito", "LICA" = "L. catesbeianus",
               "LIGR" = "L. grylio", "LISP" = "L. sphenocephalus", "NOPE" = "N. perstriatus")

month_breaks <- yday(as.Date(paste0("2001-", sprintf("%02d", 1:12), "-15")))
month_labels <- month.abb
log_breaks <- scale_y_log10(breaks = 10^seq(0, 6, 2),
                            labels = c(expression(10^0), expression(10^2), expression(10^4), expression(10^6)),
                            expand = expansion(mult = c(0, 0.08)))

# Shared layout for model-prediction curves (ribbon + line, coloured by pathogen)
curve_plot <- function(dat, x, alpha = 0.15) {
  ggplot(dat, aes({{ x }}, estimate, color = pathogen, fill = pathogen)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = alpha, color = NA) +
    geom_line(linewidth = 1.3) +
    scale_color_manual(values = path_pal, labels = path_labels, name = "Pathogen") +
    scale_fill_manual(values = path_pal, guide = "none") +
    theme_classic(base_size = 16)
}

# Marginal predictions over an environmental predictor, averaged across species x stage.
# inv = plogis for binomial models (CIs on the link scale), identity for Gaussian load models.
pred_env <- function(model, dat, predictor, x_range = range(dat[[predictor]], na.rm = TRUE),
                     inv = identity, n = 100) {
  newdat <- expand.grid(species = levels(dat$species), stage = levels(dat$stage),
                        x = seq(x_range[1], x_range[2], length.out = n)) %>%
    mutate(avg_temp_z = if (predictor == "avg_temp_z") x else 0,
           precip_z   = if (predictor == "precip_z") x else 0)
  p <- predict(model, newdata = newdat, type = "link", se.fit = TRUE, re.form = NA, allow.new.levels = TRUE)
  newdat %>%
    mutate(fit = inv(p$fit), lower = inv(p$fit - 1.96 * p$se.fit), upper = inv(p$fit + 1.96 * p$se.fit)) %>%
    group_by(x) %>%
    summarise(estimate = mean(fit), lower = mean(lower), upper = mean(upper), .groups = "drop")
}
pred_both <- function(bd_model, bd_dat, rv_model, rv_dat, predictor, ...) {
  bind_rows(Bd = pred_env(bd_model, bd_dat, predictor, ...),
            Rv = pred_env(rv_model, rv_dat, predictor, ...), .id = "pathogen") %>%
    mutate(pathogen = factor(pathogen, levels = c("Bd", "Rv")))
}

# Back-transform z-scores to observed values for plotting
temp_mean <- mean(d$avg_temp_C, na.rm = TRUE);  temp_sd <- sd(d$avg_temp_C, na.rm = TRUE)
precip_mean <- mean(d$precip_avg, na.rm = TRUE); precip_sd <- sd(d$precip_avg, na.rm = TRUE)

prev_temp <- pred_both(bd_prev_final, bd_prev_dat, rv_prev_final, rv_prev_dat, "avg_temp_z",
                       x_range = (c(5, 30) - temp_mean) / temp_sd, inv = plogis) %>%
  mutate(temp = x * temp_sd + temp_mean)
prev_precip <- pred_both(bd_prev_final, bd_prev_dat, rv_prev_final, rv_prev_dat, "precip_z",
                         x_range = (c(0, 30) - precip_mean) / precip_sd, inv = plogis) %>%
  mutate(precip = x * precip_sd + precip_mean)
load_temp   <- pred_both(bdL_final, bdL, rvL_final, rvL, "avg_temp_z") %>% mutate(temp = x * temp_sd + temp_mean)
load_precip <- pred_both(bdL_final, bdL, rvL_final, rvL, "precip_z") %>% mutate(precip = x * precip_sd + precip_mean)

# Fig. 1A: observed seasonal prevalence ----
season_prev <- d %>%
  group_by(survey_event, date, DOY, site) %>%
  summarise(n = n(), Bd = mean(bd_status, na.rm = TRUE), Rv = mean(rv_status, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(c(Bd, Rv), names_to = "pathogen", values_to = "prevalence") %>%
  mutate(pathogen = factor(pathogen, levels = c("Bd", "Rv")))

month_prev <- season_prev %>%
  mutate(month = month(date), positive = round(prevalence * n)) %>%
  group_by(pathogen, month) %>%
  summarise(n = sum(n), positive = sum(positive), prevalence = positive / n, .groups = "drop") %>%
  rowwise() %>%
  mutate(lower = binom_ci(positive, n)[1], upper = binom_ci(positive, n)[2],
         DOY = yday(as.Date(sprintf("2001-%02d-15", month)))) %>%
  ungroup()

p_prev <- ggplot() +
  geom_point(data = season_prev,
             aes(DOY, prevalence, color = pathogen, size = n),
             alpha = 0.40) +
  geom_line(data = month_prev,
            aes(DOY, prevalence, color = pathogen, group = pathogen), linewidth = 1.2) +
  geom_errorbar(data = month_prev,
                aes(DOY, ymin = lower, ymax = upper, color = pathogen), width = 5, linewidth = 0.7) +
  geom_point(data = month_prev,
             aes(DOY, prevalence, color = pathogen), size = 3.5) +
  scale_color_manual(values = path_pal,
                     name = "Pathogen",labels = path_labels) +
  scale_size_continuous(range = c(1.5, 7),
                        name = "Amphibians sampled\nper survey",
                        breaks = c(20, 40, 60),
                        labels = c("20", "40", "60+")) +
  scale_x_continuous(breaks = month_breaks,
                     labels = month_labels, limits = c(1, 365),
                     expand = expansion(mult = c(0.025, 0.025))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1),
                     expand = expansion(mult = c(0.04, 0.03))) +
  labs(x = "Month", y = "Infection prevalence") +
  guides(color = guide_legend(order = 1,
                              override.aes = list(linewidth = 1.5, size = 3, alpha = 1)), size = guide_legend(order = 2)) +
  theme_classic(base_size = 16) +
  theme(legend.position = "right")

p_prev

# Fig. 1B: seasonal cycle as odds relative to each pathogen's annual mean ----
# (logit seasonal component centered on its annual mean, exponentiated; 95% CI from 4,000 draws)
set.seed(123)
curve_sims <- MASS::mvrnorm(4000, fixef(m_full)$cond, vcov(m_full)$cond)
rel_curve <- function(beta, path) { cc <- seasonal_curves(beta)[[path]]; exp(cc - mean(cc)) }
season_rel <- map_dfr(c("Bd", "Rv"), function(path) {
  m <- apply(curve_sims, 1, rel_curve, path = path)
  tibble(DOY = doy_grid, pathogen = path, estimate = rel_curve(fixef(m_full)$cond, path),
         lower = apply(m, 1, quantile, .025), upper = apply(m, 1, quantile, .975))
}) %>% mutate(pathogen = factor(pathogen, levels = c("Bd", "Rv")))

peak_lines <- async_peaks %>% filter(grepl("peak", quantity)) %>%
  mutate(pathogen = factor(substr(quantity, 1, 2), levels = c("Bd", "Rv")))

odd_months <- c(1, 3, 5, 7, 9, 11)

p_season <- curve_plot(season_rel, DOY) +
  geom_hline(yintercept = 1, linetype = 1, linewidth = 0.6, color = "grey50") +
  geom_vline(data = peak_lines, aes(xintercept = estimate, color = pathogen),
             linetype = 3, linewidth = 0.8, show.legend = FALSE) +
  scale_x_continuous(breaks = month_breaks[odd_months], labels = month_labels[odd_months],
                     limits = c(1, 365), expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10), labels = c("0.01", "0.1", "1", "10")) +
  labs(x = "Month", y = "Relative infection odds") +
  theme(legend.position = "none")

# Fig. 1C-D: temperature and precipitation effects on infection probability ----
prev_y <- scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0, 0.03)))

p_prev_temp <- curve_plot(prev_temp, temp) + prev_y +
  scale_x_continuous(limits = c(5, 30), breaks = seq(5, 30, 5), expand = expansion(mult = c(0, 0))) +
  labs(x = "Temperature (°C)", y = "Predicted infection probability") +
  theme(legend.position = "none")

p_prev_precip <- curve_plot(prev_precip, precip) + prev_y +
  scale_x_continuous(limits = c(0, 30), breaks = c(10, 20, 30), expand = expansion(mult = c(0, 0))) +
  labs(x = "Precipitation (mm)", y = "Predicted infection probability") +
  theme(legend.position = "none")

fig1 <- p_prev / (p_season | p_prev_temp | p_prev_precip) +
  plot_layout(heights = c(1.15, 1), widths = c(1.2, 1.4), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "right")

fig1


#uncomment below to save
#ggsave("Figure1_ENV.tiff", fig1, compression = "lzw", width = 16, height = 9, units = "in", dpi = 600)

# Fig. 2A: species infection probability (HYFE Bd omitted: never Bd+) ----
prev_spp <- bind_rows(
  Bd = as.data.frame(emmeans(bd_prev_final, ~ species, type = "response")),
  Rv = as.data.frame(emmeans(rv_prev_final, ~ species, type = "response")), .id = "pathogen") %>%
  transmute(pathogen, species, estimate = prob, lower = asymp.LCL, upper = asymp.UCL) %>%
  filter(!(pathogen == "Bd" & species == "HYFE")) %>%
  mutate(pathogen = factor(pathogen, levels = c("Bd", "Rv")),
         dodge_order = factor(pathogen, levels = c("Rv", "Bd")),
         species = factor(species, levels = rev(c("ACGR", "EUQU", "HYFE", "HYGR", "LICA", "LICP", "LISP", "NOPE"))))

p_prev_species <- ggplot(prev_spp, aes(estimate, species, color = pathogen, group = dodge_order)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.45), height = 0, linewidth = 0.9) +
  geom_point(position = position_dodge(width = 0.45), shape = 16, size = 4) +
  scale_color_manual(values = path_pal, breaks = c("Bd", "Rv"), labels = path_labels, name = "Pathogen") +
  scale_y_discrete(labels = spp_names) +
  scale_x_continuous(limits = c(-0.01, 0.345), breaks = c(0, 0.1, 0.2, 0.3), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Predicted infection probability", y = NULL) +
  theme_classic(base_size = 16) +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "right")

# Fig. 2B: stage odds ratios (post vs pre) ----
stage_or <- bind_rows(
  Bd = as.data.frame(confint(pairs(emm_bd_stage), type = "response")),
  Rv = as.data.frame(confint(pairs(emm_rv_stage), type = "response")), .id = "pathogen") %>%
  transmute(pathogen = factor(pathogen, levels = c("Bd", "Rv")), odds_ratio = odds.ratio,
            lower = asymp.LCL, upper = asymp.UCL)

p_prev_stage_or <- ggplot(stage_or, aes(odds_ratio, pathogen, color = pathogen)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.7, color = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.9) +
  geom_point(shape = 16, size = 4) +
  scale_color_manual(values = path_pal, labels = path_labels, name = "Pathogen") +
  scale_y_discrete(labels = path_labels, limits = rev) +
  scale_x_continuous(breaks = 0:3, limits = c(0.7, 3)) +
  labs(x = "Infection odds ratio (Post vs. Pre-metamorphic)", y = NULL) +
  theme_classic(base_size = 16)

# Fig. 2C-D: loads by species and stage (geometric mean +/- 95% CI) ----
load_plot_dat <- bind_rows(
  bdL %>% transmute(species, metam_stage = stage, pathogen = "Bd", load = bd_its + 1),
  rvL %>% transmute(species, metam_stage = stage, pathogen = "Rv", load = rv_ge + 1)) %>%
  mutate(pathogen = factor(pathogen, levels = c("Bd", "Rv")),
         species = factor(species, levels = c("ACGR", "HYFE", "HYGR", "EUQU", "LICP", "LICA", "LISP", "NOPE")))

geom_summary <- function(dat, ...) {
  dat %>% filter(load > 0) %>% mutate(log_load = log10(load)) %>%
    group_by(...) %>%
    summarise(n = n(), mean_log = mean(log_load), se_log = sd(log_load) / sqrt(n),
              tcrit = qt(0.975, df = pmax(n - 1, 1)), geometric_mean = 10^mean_log,
              lower = 10^(mean_log - tcrit * se_log), upper = 10^(mean_log + tcrit * se_log), .groups = "drop")
}

load_species_summary <- geom_summary(load_plot_dat, species, pathogen)

pd <- position_dodge(width = 0.65)

p_load_species <- ggplot(load_plot_dat, aes(species, load, color = pathogen)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.65, seed = 123),
             alpha = 0.25, size = 2.5) +
  geom_errorbar(data = load_species_summary, aes(species, ymin = lower, ymax = upper, group = pathogen),
                position = pd, width = 0, linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(data = load_species_summary, aes(species, geometric_mean, group = pathogen),
             position = pd, size = 3.6, shape = 16, color = "black", inherit.aes = FALSE) +
  scale_color_manual(values = path_pal, name = "Pathogen", labels = path_labels) +
  scale_x_discrete(labels = spp_names) +
  log_breaks + coord_cartesian(ylim = c(0.5, NA)) +
  labs(x = "Host species", y = "Log10 pathogen load") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11), legend.position = "right")

stage_levels <- c("Bd_pre", "Bd_post", "Rv_pre", "Rv_post")
load_plot_stage <- load_plot_dat %>% filter(!is.na(metam_stage)) %>%
  mutate(stage_group = factor(paste(pathogen, metam_stage, sep = "_"), levels = stage_levels))
load_stage_summary <- geom_summary(load_plot_stage, metam_stage, pathogen) %>%
  mutate(stage_group = factor(paste(pathogen, metam_stage, sep = "_"), levels = stage_levels))

p_load_stage <- ggplot(load_plot_stage, aes(stage_group, load, color = pathogen)) +
  ggdist::stat_halfeye(aes(fill = pathogen), adjust = 0.7, width = 0.55, justification = -0.35,
                       slab_alpha = 0.22, point_interval = NULL, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 123), alpha = 0.25, size = 2.5) +
  geom_errorbar(data = load_stage_summary, aes(stage_group, ymin = lower, ymax = upper),
                color = "black", width = 0, linewidth = 1, inherit.aes = FALSE) +
  geom_point(data = load_stage_summary, aes(stage_group, geometric_mean), color = "black", size = 3.6, inherit.aes = FALSE) +
  scale_color_manual(values = path_pal, name = "Pathogen", labels = path_labels) +
  scale_fill_manual(values = path_pal, guide = "none") +
  scale_x_discrete(labels = c("Bd_pre" = "Pre", "Bd_post" = "Post", "Rv_pre" = "Pre", "Rv_post" = "Post")) +
  log_breaks + coord_cartesian(ylim = c(0.5, NA)) +
  labs(x = "Metamorphic stage", y = "Log10 pathogen load") +
  theme_classic(base_size = 16) +
  theme(legend.position = "right", axis.ticks.x = element_blank())

# Panel A: keeps the only legend, drawn as solid squares
p_prev_species <- p_prev_species +
  guides(color = guide_legend(override.aes = list(linetype = 0, shape = 16, size = 5, alpha = 1)))

# Panels B-D: drop their copies
p_prev_stage_or <- p_prev_stage_or + guides(color = "none")
p_load_species  <- p_load_species  + guides(color = "none")
p_load_stage    <- p_load_stage    + guides(color = "none")

fig2 <- (p_prev_species | p_prev_stage_or) / (p_load_species | p_load_stage) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right", legend.text = element_text(hjust = 0))
fig2


ggsave("Figure2_HOST.tiff", fig2, compression = "lzw", width = 15, height = 10, units = "in", dpi = 600)

# Fig. S1: temperature and precipitation effects on load ----
load_y <- scale_y_continuous(breaks = 0:6, labels = scales::label_math(10^.x), expand = expansion(mult = c(0, 0.05)))

p_load_temp <- curve_plot(load_temp, temp, alpha = 0.20) + load_y + coord_cartesian(ylim = c(0, 6)) +
  labs(x = "Temperature (°C)", y = "Predicted infection load (log10)") +
  theme(legend.position = "none")
p_load_precip <- curve_plot(load_precip, precip, alpha = 0.20) + load_y + coord_cartesian(ylim = c(0, 6)) +
  labs(x = "Precipitation (mm)", y = "Predicted infection load (log10)") +
  theme(legend.position = "right")

figS1 <- p_load_temp + p_load_precip + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")
figS1
#ggsave("Figure_S1_load.tiff", figS1, compression = "lzw", width = 11.5, height = 4.5, units = "in", dpi = 600)

# Fig. S2: coinfection by stage (A) and permutation results (B) ----
coinf_stage <- d %>%
  filter(species %in% c("ACGR", "NOPE"), !is.na(bd_status), !is.na(rv_status)) %>%
  mutate(coinfected = bd_status == 1 & rv_status == 1) %>%
  group_by(species, stage) %>%
  summarise(n = n(), Bd = mean(bd_status), Rv = mean(rv_status), Coinfected = mean(coinfected), .groups = "drop") %>%
  pivot_longer(c(Bd, Rv, Coinfected), names_to = "infection", values_to = "prevalence") %>%
  mutate(species = factor(species, levels = c("ACGR", "NOPE")),
         stage = factor(stage, levels = c("pre", "post")),
         infection = factor(infection, levels = c("Bd", "Rv", "Coinfected")))

p_coinf_stage <- ggplot(coinf_stage, aes(stage, prevalence, color = infection, group = infection)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3.5) +
  facet_wrap(~species, labeller = labeller(species = c("ACGR" = "Acris gryllus", "NOPE" = "Notophthalmus perstriatus"))) +
  scale_color_manual(values = c(path_pal, "Coinfected" = "#44AA99"), name = "Infection status",
                     labels = c(path_labels, "Coinfected" = "Coinfected")) +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1), labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Metamorphic stage", y = "Prevalence", color = NULL) +
  theme_classic(base_size = 16) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "italic"), legend.position = "top")

p_coinf_perm <- coinf_perm_results %>%
  mutate(species = factor(species, levels = c("ACGR", "NOPE"), labels = c("Acris gryllus", "Notophthalmus perstriatus"))) %>%
  ggplot(aes(species)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, linewidth = 0.8, alpha = 0.5) +
  geom_point(aes(y = expected, shape = "Expected"), size = 3.5, alpha = 0.5) +
  geom_point(aes(y = observed, shape = "Observed"), size = 4, color = "darkred") +
  scale_shape_manual(values = c("Expected" = 16, "Observed" = 17), name = NULL) +
  scale_y_continuous(limits = c(15, 40), breaks = seq(15, 40, 5)) +
  labs(x = NULL, y = "Number of coinfections") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(face = "italic"), legend.position = "top")

figS2 <- (p_coinf_stage | p_coinf_perm) + plot_annotation(tag_levels = "A") + plot_layout(widths = c(1.3, 1.3))
figS2
#ggsave("Figure_S2_coinfection.tiff", figS2, width = 12, height = 5, units = "in", compression = "lzw", dpi = 600)
