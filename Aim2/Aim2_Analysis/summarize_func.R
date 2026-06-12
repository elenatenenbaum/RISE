# ------------------------------------------------------------
# make_flow_table()
#
# Builds a participant flow table showing attrition at each
# stage of data processing, broken down by Age group.
#
# Stages:
#   1. Enrolled        — all consented participants (from Part_Age)
#   2. iCatcher_Ran    — participants with iCatcher frame counts
#                        (% of enrolled)
#   3. iCatcher_Good   — passed iCatcher quality check
#                        (% of iCatcher ran)
#   4. Fixation        — excluded at fixation trial stage
#                        (% of iCatcher good)
#   5. Usable          — final analytic sample
#
# Inputs:
#   results    : output of process_icatcher(); uses frame_counts
#   enrollment : list with metrics (final sample) and
#                metrics_1 (post-iCatcher, pre-fixation)
#   Part_Age   : consented analytic sample with Age column
#
# Returns:
#   A data frame with one row per Age group and columns:
#   Age, Enrolled, iCatcher_Ran, iCatcher_Good, Fixation, Usable
# ------------------------------------------------------------

make_flow_table <- function(results, enrollment, Part_Age) {
  
  enrolled <- Part_Age %>%
    group_by(Age) %>%
    summarise(n_enrolled = n(), .groups = "drop")
  
  icatcher_ran <- results$frame_counts %>%
    merge(Part_Age, all.x = TRUE) %>%
    distinct(response_uuid, .keep_all = TRUE) %>%
    group_by(Age) %>%
    summarise(n_icatcher = n(), .groups = "drop")
  
  icatcher_good <- enrollment$metrics_1 %>%
    distinct(response_uuid, .keep_all = TRUE) %>%
    group_by(Age) %>%
    summarise(n_icatcher_good = n(), .groups = "drop")
  
  fixation_excluded <- enrollment$metrics_1 %>%
    distinct(response_uuid, .keep_all = TRUE) %>%
    filter(!response_uuid %in% enrollment$metrics$response_uuid) %>%
    group_by(Age) %>%
    summarise(n_fixation = n(), .groups = "drop")
  
  usable <- enrollment$metrics %>%
    distinct(response_uuid, .keep_all = TRUE) %>%
    group_by(Age) %>%
    summarise(n_usable = n(), .groups = "drop")
  
  enrolled %>%
    left_join(icatcher_ran,       by = "Age") %>%
    left_join(icatcher_good,      by = "Age") %>%
    left_join(fixation_excluded,  by = "Age") %>%
    left_join(usable,             by = "Age") %>%
    mutate(
      Enrolled      = as.character(n_enrolled),
      iCatcher_Ran  = sprintf("%d (%.1f%%)", n_icatcher,      (n_icatcher      / n_enrolled)   * 100),
      iCatcher_Good = sprintf("%d (%.1f%%)", n_icatcher_good, (n_icatcher_good / n_icatcher)   * 100),
      Fixation      = sprintf("%d (%.1f%%)", n_fixation,      (n_fixation      / n_icatcher_good) * 100),
      Usable        = n_usable
    ) %>%
    dplyr::select(Age, Enrolled, iCatcher_Ran, iCatcher_Good, Fixation, Usable)
}

# ------------------------------------------------------------
# summarise_quality_v2()
#
# Computes per-participant video quality metrics and returns
# a gtsummary table stratified by Age group.
#
# Quality metrics (averaged across trials per participant):
#   - noface_p  : % frames with no face detected
#   - away_p    : % frames coded as away looks
#   - lowconf_p : % left/right frames below 85% confidence
#
# All percentages are calculated out of frames with iCatcher+
# output only (n_icatcher = total frames - missing frames).
# Trials with no iCatcher output are counted as missing trials.
#
# Inputs:
#   results  : output of process_icatcher(); uses frame_counts
#              and quality_flags
#   Part_Age : participant metadata with Age, cohort, Gender
#
# Returns:
#   A gtsummary tbl_summary object 
# ------------------------------------------------------------

summarise_quality_v2 <- function(results, Part_Age) {
  
  base_data <- results$frame_counts %>%
    merge(results$quality_flags %>% dplyr::select(-Age),
          by  = c("response_uuid", "Trial"),
          all = TRUE) %>%
    mutate(noicatcherdata = tidyr::replace_na(noicatcherdata, 0),
           n_icatcher     = n - noicatcherdata)
  
  per_participant <- base_data %>%
    mutate(
      away    = tidyr::replace_na(away,    0),
      noface  = tidyr::replace_na(noface,  0),
      lowconf = tidyr::replace_na(lowconf, 0),
      noface_p  = ifelse(n_icatcher > 0, (noface  / n_icatcher) * 100, NA),
      away_p    = ifelse(n_icatcher > 0, (away    / n_icatcher) * 100, NA),
      lowconf_p = ifelse(n_icatcher > 0, (lowconf / n_icatcher) * 100, NA)
    ) %>%
    group_by(response_uuid) %>%
    dplyr::summarise(
      max_trial = max(n()),
      missing_trial = sum(n_icatcher == 0),
      has_missing   = missing_trial > 0,
      noface_p      = mean(noface_p,  na.rm = TRUE),
      away_p        = mean(away_p,    na.rm = TRUE),
      lowconf_p     = mean(lowconf_p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(c(noface_p, away_p, lowconf_p), ~ ifelse(is.nan(.), NA, .))) %>%
    merge(Part_Age, all.x = TRUE)
  
  # Table 1: quality metrics for all participants
  t1 <- per_participant %>%
    dplyr::select(-response_uuid, -cohort, -Gender, -has_missing, -missing_trial) %>%
    tbl_summary(
      by        = Age,
      type      = list(noface_p:lowconf_p ~ "continuous"),
      statistic = list(noface_p:lowconf_p ~ "{median} ({p25}, {p75})"),
      digits    = all_continuous() ~ 1,
      label     = list(
        noface_p  ~ "No face detected (%)",
        away_p    ~ "Away looks (%)",
        lowconf_p ~ "Low confidence left & right frames (%)"
      )
    ) %>%
    modify_caption(
      "Video quality metrics by age group. Percentages reflect the mean across trials 
       with iCatcher+ output only, calculated out of frames processed by iCatcher+. 
       Trials where iCatcher+ produced no output are counted as missing trials."
    )
  
  t1
}

# ------------------------------------------------------------
# summarise_trials()
#
# Summarises the number of usable test trials per participant
# in the final analytic sample, stratified by Age group.
#
# Returns median (IQR) and min/max trial counts as a
# gtsummary table.
#
# Inputs:
#   enrollment : list with metrics (final analytic sample)
#                containing response_uuid, Age, and one row
#                per usable trial
#
# Returns:
#   A gtsummary tbl_summary object
# ------------------------------------------------------------

summarise_trials <- function(enrollment) {
  enrollment$metrics %>%
    dplyr::group_by(response_uuid) %>%
    dplyr::mutate(trial_num = n()) %>%
    distinct(response_uuid, .keep_all = TRUE) %>%
    ungroup() %>%
    dplyr::select(trial_num, Age) %>%
    tbl_summary(
      by   = Age,
      type = list(trial_num ~ "continuous2"),
      statistic = list(
        all_continuous()  ~ c("{median} ({p25}, {p75})", "{min}, {max}"),
        all_categorical() ~ "{n} / {N} ({p}%)"
      ),
      digits       = all_continuous() ~ 2,
      missing_text = "(Missing)"
    ) %>% modify_caption("**Median (IQR) number of test trials with usable data (final analytic)**")
}
