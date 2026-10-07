# =============================================================================
# PD Warrior Brain Task- LMM Models
#
#   Primary      log(outcome) ~ exercise + (1 + exercise | participant)
#   Repetition   log(outcome) ~ exercise * repetition
#                               + (1 + exercise * repetition | participant)
#   RML          log(outcome) ~ exercise * RML_difference
#                               + (1 + exercise | participant)
#   Adjusted     log(outcome) ~ exercise + age + disease duration
#                               + RML_difference + (1 + exercise | participant)
#
# Each fitted separately for speed and dwell time in each hand.
# All by REML in lme4, Satterthwaite degrees of freedom via lmerTest.
# Confidence intervals use the t distribution on the Satterthwaite df.
#
# Leave-one-out analysis performed for the primary models, 
# Relative changes in the exercise estimate on adding predictors calculated 
# singularity checked for every model.

# =============================================================================

library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)

base_dir <- paste0(
  "C:/Users/lz23950/OneDrive - University of Bristol/",
  "brain_task_plots/actual_code/"
)

out_dir <- "."


# =============================================================================
# 1. FILE AND COLUMN MAPPING
# =============================================================================

outcome_spec <- list(
  speed = list(
    file       = paste0(base_dir, "velocity_more_standardized_as_continous_FOR_R.csv"),
    id         = "IDrep2",
    exercise   = "exerep2",
    repetition = "recrep2",
    value      = "allV",
    label      = "Speed"
  ),
  dwell = list(
    file       = paste0(base_dir, "dwell_more_standardized_as_continous_FOR_R.csv"),
    id         = "IDrep1",
    exercise   = "exerep1",
    repetition = "recrep1",
    value      = "allDwell",
    label      = "Dwell time"
  )
)

HAND_COL      <- "RecordingWithMoreAffectedHand"
CONTROL_LEVEL <- 0

# Age, disease duration and relative medication level are participant-level
# constants, so they are read once from the speed file and joined to both
# datasets by participant ID.
PERSON_FILE <- "speed"
AGE_COL     <- "AGE_SD"
DD_COL      <- "DD_SD"
RML_COL     <- "RelativePlasmaLevodopa_standardized"

raw_data <- lapply(outcome_spec, function(s) {
  d <- read.csv(s$file, stringsAsFactors = FALSE)
  need <- c(s$id, s$exercise, s$repetition, s$value, HAND_COL)
  miss <- setdiff(need, names(d))
  if (length(miss) > 0) {
    cat("\nMissing in ", basename(s$file), ": ",
        paste(miss, collapse = ", "), "\nAvailable columns:\n", sep = "")
    print(names(d))
    stop("Fix the column mapping in Section 1.")
  }
  d
})

need_person <- c(AGE_COL, DD_COL, RML_COL)
miss_person <- setdiff(need_person, names(raw_data[[PERSON_FILE]]))
if (length(miss_person) > 0) {
  cat("\nMissing participant-level columns in ",
      basename(outcome_spec[[PERSON_FILE]]$file), ": ",
      paste(miss_person, collapse = ", "), "\nAvailable columns:\n", sep = "")
  print(names(raw_data[[PERSON_FILE]]))
  stop("Fix AGE_COL / DD_COL / RML_COL in Section 1.")
}


# =============================================================================
# 2. PARTICIPANT-LEVEL RML DIFFERENCE
# =============================================================================
#
# One RML value per participant per condition; the difference is
# exercise minus control, z-scored.

sv <- outcome_spec$speed

person_long <- raw_data[[PERSON_FILE]] %>%
  transmute(
    ID       = as.character(.data[[sv$id]]),
    Exercise = ifelse(.data[[sv$exercise]] == CONTROL_LEVEL, "Control", "Exercise"),
    RML      = as.numeric(.data[[RML_COL]]),
    Age_z    = as.numeric(.data[[AGE_COL]]),
    DD_z     = as.numeric(.data[[DD_COL]])
  )

person_fixed <- person_long %>%
  group_by(ID) %>%
  summarise(Age_z = first(Age_z), DD_z = first(DD_z), .groups = "drop")

person <- person_long %>%
  group_by(ID, Exercise) %>%
  summarise(RML = first(RML), .groups = "drop") %>%
  pivot_wider(names_from = Exercise, values_from = RML, names_prefix = "RML_") %>%
  mutate(
    RML_diff   = RML_Exercise - RML_Control,
    RML_diff_z = as.numeric(scale(RML_diff))
  ) %>%
  left_join(person_fixed, by = "ID")

cat("\n============================================================\n")
cat("RML SESSION DIFFERENCE (participant level)\n")
cat("============================================================\n")
print(as.data.frame(person), row.names = FALSE, digits = 4)
cat("\nSD used for z-scoring:", round(sd(person$RML_diff), 4), "\n")


# =============================================================================
# 3. DATA PREPARATION
# =============================================================================

prepare <- function(outcome_name, hand_value) {

  s <- outcome_spec[[outcome_name]]

  d <- raw_data[[outcome_name]] %>%
    filter(.data[[HAND_COL]] == hand_value) %>%
    transmute(
      ID         = as.character(.data[[s$id]]),
      Exercise   = factor(.data[[s$exercise]],
                          levels = c(CONTROL_LEVEL, 1 - CONTROL_LEVEL),
                          labels = c("Control", "Exercise")),
      Repetition = factor(.data[[s$repetition]]),
      y          = log(as.numeric(.data[[s$value]]))
    ) %>%
    filter(is.finite(y)) %>%
    left_join(person[, c("ID", "RML_diff_z", "Age_z", "DD_z")], by = "ID") %>%
    filter(complete.cases(.)) %>%
    mutate(ID = factor(ID), Repetition = droplevels(Repetition))

  contrasts(d$Exercise) <- matrix(
    c(-0.5, 0.5), ncol = 1, dimnames = list(NULL, "ExVsCon"))
  contrasts(d$Repetition) <- matrix(
    c(-0.5, 0.5), ncol = 1, dimnames = list(NULL, "R2vsR1"))
  d
}

combos <- expand.grid(outcome = names(outcome_spec), hand = c(1, 0),
                      stringsAsFactors = FALSE)
combos$hand_label <- ifelse(combos$hand == 1, "More affected", "Less affected")
combos$label <- paste(sapply(combos$outcome, function(o) outcome_spec[[o]]$label),
                      combos$hand_label, sep = ", ")
combos <- combos[order(combos$outcome, -combos$hand), ]

# =============================================================================
# 4. EXTRACTION HELPER
# =============================================================================

tidy_fit <- function(fit, outcome_label, model_label) {
  cf <- summary(fit)$coefficients
  data.frame(
    outcome  = outcome_label,
    model    = model_label,
    term     = rownames(cf),
    df       = cf[, "df"],
    beta     = cf[, "Estimate"],
    SE       = cf[, "Std. Error"],
    t        = cf[, "t value"],
    p        = cf[, "Pr(>|t|)"],
    exp_beta = exp(cf[, "Estimate"]),
    CI_low   = exp(cf[, "Estimate"] - qt(0.975, cf[, "df"]) * cf[, "Std. Error"]),
    CI_high  = exp(cf[, "Estimate"] + qt(0.975, cf[, "df"]) * cf[, "Std. Error"]),
    singular = isSingular(fit, tol = 1e-5),
    row.names = NULL
  )
}

pretty <- function(x) {
  x <- gsub("ExerciseExVsCon",   "Exercise", x)
  x <- gsub("RepetitionR2vsR1",  "Repetition", x)
  x <- gsub("Age_z",             "Age", x)
  x <- gsub("DD_z",              "Disease duration", x)
  x <- gsub("RML_diff_z",        "RML difference", x)
  gsub(":", " x ", x)
}

fmt <- function(df) {
  df %>% mutate(
    df       = round(df, 1),
    beta     = round(beta, 4),
    SE       = round(SE, 4),
    t        = round(t, 2),
    p        = ifelse(p < .001, "< .001", sub("^0", "", sprintf("%.3f", p))),
    exp_beta = round(exp_beta, 3),
    CI       = sprintf("[%.3f, %.3f]", CI_low, CI_high),
    term     = pretty(term)
  ) %>%
    select(outcome, model, term, df, beta, SE, t, p, exp_beta, CI)
}


# =============================================================================
# 5. FIT THE FOUR MODEL TYPES
# =============================================================================

formulas <- list(
  Primary    = "y ~ Exercise + (1 + Exercise | ID)",
  Repetition = "y ~ Exercise * Repetition + (1 + Exercise * Repetition | ID)",
  RML        = "y ~ Exercise * RML_diff_z + (1 + Exercise | ID)",
  Adjusted   = "y ~ Exercise + Age_z + DD_z + RML_diff_z + (1 + Exercise | ID)"
)

all_estimates <- list()
exercise_beta <- list()

for (i in seq_len(nrow(combos))) {

  d <- prepare(combos$outcome[i], combos$hand[i])

  for (mod in names(formulas)) {
    fit <- suppressWarnings(
      lmerTest::lmer(as.formula(formulas[[mod]]), data = d, REML = TRUE))
    all_estimates[[length(all_estimates) + 1]] <-
      tidy_fit(fit, combos$label[i], mod)
    exercise_beta[[length(exercise_beta) + 1]] <- data.frame(
      outcome = combos$label[i], model = mod,
      beta = unname(fixef(fit)["ExerciseExVsCon"]))
  }
}

estimates     <- bind_rows(all_estimates)
exercise_beta <- bind_rows(exercise_beta)


# ---- Table 2: primary models ------------------------------------------------

cat("\n\n============================================================\n")
cat("TABLE 2 — PRIMARY MODELS\n")
cat("============================================================\n")
print(as.data.frame(
  fmt(estimates %>% filter(model == "Primary", term != "(Intercept)"))),
  row.names = FALSE)


# ---- Supplementary Table 1: adjusted models ---------------------------------

cat("\n\n============================================================\n")
cat("SUPPLEMENTARY TABLE 1 — ADJUSTED MODELS\n")
cat("============================================================\n")
print(as.data.frame(
  fmt(estimates %>% filter(model == "Adjusted", term != "(Intercept)"))),
  row.names = FALSE)


# ---- Supplementary Table 2: interaction models ------------------------------

cat("\n\n============================================================\n")
cat("SUPPLEMENTARY TABLE 2 — INTERACTION MODELS (REML)\n")
cat("============================================================\n")
print(as.data.frame(
  fmt(estimates %>%
        filter(model %in% c("Repetition", "RML"), term != "(Intercept)"))),
  row.names = FALSE)

cat("\nThese are REML fits, as the Methods state. The values currently in\n")
cat("Supplementary Table 2 came from ML fits; replace them with these.\n")


# =============================================================================
# 6. LEAVE-ONE-OUT ON THE PRIMARY MODELS
# =============================================================================

loo <- bind_rows(lapply(seq_len(nrow(combos)), function(i) {

  d <- prepare(combos$outcome[i], combos$hand[i])
  ids <- levels(d$ID)

  vals <- sapply(ids, function(drop) {
    dd <- droplevels(d[d$ID != drop, ])
    f  <- suppressWarnings(
      lmerTest::lmer(y ~ Exercise + (1 + Exercise | ID), data = dd, REML = TRUE))
    unname(fixef(f)["ExerciseExVsCon"])
  })

  full <- unname(fixef(suppressWarnings(
    lmerTest::lmer(y ~ Exercise + (1 + Exercise | ID),
                   data = d, REML = TRUE)))["ExerciseExVsCon"])

  data.frame(
    outcome   = combos$label[i],
    full_beta = full,
    full_pct  = (exp(full) - 1) * 100,
    loo_min_pct = (exp(min(vals)) - 1) * 100,
    loo_max_pct = (exp(max(vals)) - 1) * 100,
    n_fits    = length(vals)
  )
}))

cat("\n\n============================================================\n")
cat("LEAVE-ONE-OUT RANGES — PRIMARY MODELS\n")
cat("============================================================\n")
print(as.data.frame(loo), row.names = FALSE, digits = 4)


# =============================================================================
# 7. RELATIVE CHANGE IN THE EXERCISE ESTIMATE
# =============================================================================


pct_change <- exercise_beta %>%
  pivot_wider(names_from = model, values_from = beta) %>%
  mutate(
    pct_change_Repetition = 100 * abs(Repetition - Primary) / abs(Primary),
    pct_change_RML        = 100 * abs(RML - Primary)        / abs(Primary),
    pct_change_Adjusted   = 100 * abs(Adjusted - Primary)   / abs(Primary)
  )

cat("\n\n============================================================\n")
cat("RELATIVE CHANGE IN THE EXERCISE ESTIMATE (%)\n")
cat("============================================================\n")
print(as.data.frame(
  pct_change[, c("outcome", "Primary", "Repetition", "RML", "Adjusted",
                 "pct_change_Repetition", "pct_change_RML",
                 "pct_change_Adjusted")]),
  row.names = FALSE, digits = 4)

cat("\nMaximum change on adding repetition :",
    round(max(pct_change$pct_change_Repetition), 2), "%\n")
cat("Maximum change on adding RML         :",
    round(max(pct_change$pct_change_RML), 2), "%\n")
cat("Maximum change on adding covariates  :",
    round(max(pct_change$pct_change_Adjusted), 2), "%\n")
cat("Maximum across all three             :",
    round(max(c(pct_change$pct_change_Repetition,
                pct_change$pct_change_RML,
                pct_change$pct_change_Adjusted)), 2), "%\n")


# =============================================================================
# 8. CONVERGENCE AND SINGULARITY
# =============================================================================

cat("\n\n============================================================\n")
cat("SINGULARITY CHECK\n")
cat("============================================================\n")
sing <- estimates %>%
  group_by(outcome, model) %>%
  summarise(singular = any(singular), .groups = "drop")
print(as.data.frame(sing), row.names = FALSE)
cat("\nAny TRUE invalidates the statement that all models converged without\n")
cat("singularity.\n")


# =============================================================================
# 9. SAVE
# =============================================================================

write.csv(fmt(estimates), file.path(out_dir, "all_model_estimates.csv"),
          row.names = FALSE)
write.csv(loo,        file.path(out_dir, "leave_one_out_ranges.csv"), row.names = FALSE)
write.csv(pct_change, file.path(out_dir, "relative_change.csv"),      row.names = FALSE)

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
