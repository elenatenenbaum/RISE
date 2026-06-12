# ==============================================================================
# TABLE 3 & TABLE 4
# Drop this chunk after all task sections are complete.
# ==============================================================================

# ------------------------------------------------------------------------------
# TABLE 3 HELPER
# Computes one row of Table 3 (data quality metrics) for a given task.
#
# For dual-screen tasks (process_icatcher):   quality_flags has away, noface, lowconf, noicatcherdata
# For single-screen AV (process_icatcher_single): quality_flags has noface, lowconf, noicatcherdata (no away)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# TABLE 4 HELPER
# Computes one row of Table 4 (participant retention) for a given task.
#
# Side bias N = participants in metrics_1 who are NOT in metrics.
# For tasks with remove_fixations = FALSE (Prediction, AV), this will be 0
# by definition since no fixation filtering was applied — which is correct.
#
# trial_col: name of the trial ID column in enrollment$metrics.
#   Most tasks: "Trial"
#   VPC and HH (collapse_vpc_hh = TRUE): "Trial Total"
# ------------------------------------------------------------------------------


# ==============================================================================
# BUILD TABLE 4
# N valid iCatcher per task comes from the icatcher_ran step in make_flow_table
# (frame_counts distinct by response_uuid). We replicate that here.
# ==============================================================================

get_n_icatcher <- function(results, Part_Age) {
  results$frame_counts %>%
    merge(Part_Age, all.x = TRUE) %>%
    distinct(response_uuid) %>%
    nrow()
}


# ------------------------------------------------------------------------------
# HELPER: compute per-age-group string "N (pct%)" for any count vector
# ------------------------------------------------------------------------------

age_str <- function(ns, pcts) {
  paste(mapply(function(n, p) paste0(n, " (", p, "%)"), ns, pcts), collapse = " / ")
}

# ------------------------------------------------------------------------------
# TABLE 3 HELPER AIM 2
# ------------------------------------------------------------------------------

compute_table3_row <- function(results, icatcher_df, task_label,
                               n_trials_qc, Part_Age,
                               show_away   = TRUE,
                               quality_tbl = NULL) {
  
  age_levels  <- sort(unique(Part_Age$Age))
  age_suffix  <- gsub("[^0-9]", "", age_levels)
  by_levels   <- age_levels  # just alias, no separate sort()
  age_counts  <- Part_Age %>%
    group_by(Age) %>%
    summarise(n_consented = n(), .groups = "drop") %>%
    arrange(Age)
  
  n_consented <- nrow(Part_Age)
  
  # --- Valid iCatcher N overall and by age ---
  icatcher_ages <- icatcher_df %>%
    merge(Part_Age %>% dplyr::select(response_uuid, Age),
          by = "response_uuid", all.x = TRUE) %>%
    distinct(response_uuid, Age)
  
  valid_by_age <- icatcher_ages %>%
    group_by(Age) %>%
    summarise(n_valid = n(), .groups = "drop") %>%
    merge(age_counts, by = "Age") %>%
    mutate(pct_valid = sprintf("%.1f", n_valid / n_consented * 100)) %>%
    arrange(Age)
  
  n_valid   <- length(unique(icatcher_df$response_uuid))
  pct_valid <- sprintf("%.1f", n_valid / n_consented * 100)
  
  # --- N with 0% task completion overall and by age ---
  base_data <- results$frame_counts %>%
    merge(results$quality_flags %>% dplyr::select(-Age),
          by  = c("response_uuid", "Trial"),
          all = TRUE) %>%
    mutate(noicatcherdata = tidyr::replace_na(noicatcherdata, 0),
           n_icatcher     = n - noicatcherdata)
  
  zero_ids <- base_data %>%
    group_by(response_uuid) %>%
    summarise(all_zero = all(n_icatcher == 0), .groups = "drop") %>%
    filter(all_zero) %>%
    pull(response_uuid)
  
  zero_by_age <- data.frame(response_uuid = zero_ids) %>%
    merge(Part_Age %>% dplyr::select(response_uuid, Age),
          by = "response_uuid", all.x = TRUE) %>%
    group_by(Age) %>%
    summarise(n_zero = n(), .groups = "drop") %>%
    merge(data.frame(Age = age_levels), all.y = TRUE) %>%
    mutate(n_zero = tidyr::replace_na(n_zero, 0)) %>%
    merge(valid_by_age %>% dplyr::select(Age, n_valid), by = "Age") %>%
    mutate(pct_zero = sprintf("%.1f", n_zero / n_valid * 100)) %>%
    arrange(Age)
  
  zero_completion <- length(zero_ids)
  pct_zero        <- sprintf("%.1f", zero_completion / n_valid * 100)
  
  # --- Quality table ---
  if (is.null(quality_tbl)) {
    quality_tbl <- summarise_quality_v2(results, Part_Age)
  }
  
  tbl_body  <- quality_tbl$table_body
  stat_cols <- grep("^stat_", colnames(tbl_body), value = TRUE)
  by_levels <- sort(age_levels)
  
  extract_stat <- function(var_name) {
    row <- tbl_body %>%
      filter(variable == var_name, row_type == "label")
    if (nrow(row) == 0) return(rep(NA_character_, length(stat_cols)))
    vals <- row %>% dplyr::select(all_of(stat_cols)) %>% unlist()
    trimws(vals)
  }
  
  noface_vals  <- extract_stat("noface_p")
  away_vals    <- if (show_away) extract_stat("away_p") else rep("---", length(by_levels))
  lowconf_vals <- extract_stat("lowconf_p")
  
  # --- Build output row with split columns ---
  row <- data.frame(
    Task            = task_label,
    N_trials_QC     = n_trials_qc,
    stringsAsFactors = FALSE
  )
  
  # Add per-age columns dynamically
  for (i in seq_along(by_levels)) {
    sfx <- age_suffix[i]
    vba <- valid_by_age %>% filter(Age == by_levels[i])
    zba <- zero_by_age  %>% filter(Age == by_levels[i])
    
    row[[paste0("N_valid_", sfx)]] <- paste0(vba$n_valid, " (", vba$pct_valid, "%)")
    row[[paste0("N_zero_",  sfx)]] <- paste0(zba$n_zero,  " (", zba$pct_zero,  "%)")
  }
  
  for (i in seq_along(by_levels)) {
    sfx <- age_suffix[i]
    row[[paste0("Med_pct_noface_",  sfx)]] <- noface_vals[i]
    row[[paste0("Med_pct_away_",    sfx)]] <- away_vals[i]
    row[[paste0("Med_pct_lowconf_", sfx)]] <- lowconf_vals[i]
  }
  
  row
}

# ------------------------------------------------------------------------------
# TABLE 4 HELPER
# ------------------------------------------------------------------------------

compute_table4_row <- function(enrollment, task_label,
                               n_analytic_trials,
                               Part_Age,
                               trial_col = "Trial") {
  
  age_levels <- sort(unique(Part_Age$Age))
  
  n_usable_1trial <- length(unique(enrollment$metrics_1$response_uuid))
  
  ids_post_usable   <- unique(enrollment$metrics_1$response_uuid)
  ids_post_fixation <- unique(enrollment$metrics$response_uuid)
  n_sidebias        <- length(setdiff(ids_post_usable, ids_post_fixation))
  pct_sidebias      <- sprintf("%.1f", n_sidebias / n_usable_1trial * 100)
  
  final_n <- length(unique(enrollment$metrics$response_uuid))
  
  if (!trial_col %in% colnames(enrollment$metrics)) {
    warning(paste("trial_col =", trial_col, "not found for", task_label,
                  "— defaulting to 'Trial'"))
    trial_col <- "Trial"
  }
  
  # Trials per person with age attached
  trials_per_person <- enrollment$metrics %>%
    dplyr::select(-any_of("Age")) %>%
    merge(Part_Age %>% dplyr::select(response_uuid, Age),
          by = "response_uuid", all.x = TRUE) %>%
    group_by(response_uuid, Age) %>%
    summarise(n_trials = n_distinct(.data[[trial_col]]), .groups = "drop")
  
  # Overall median IQR — rounded to 1 dp
  t_med <- round(median(trials_per_person$n_trials), 1)
  t_q1  <- round(quantile(trials_per_person$n_trials, 0.25), 1)
  t_q3  <- round(quantile(trials_per_person$n_trials, 0.75), 1)
  
  # By-age median IQR — rounded to 1 dp
  trials_by_age <- trials_per_person %>%
    group_by(Age) %>%
    summarise(
      t_med = round(median(n_trials), 1),
      t_q1  = round(quantile(n_trials, 0.25), 1),
      t_q3  = round(quantile(n_trials, 0.75), 1),
      n     = n(),
      .groups = "drop"
    ) %>%
    arrange(Age)
  
  # By-age usable and sidebias counts
  usable_by_age <- enrollment$metrics_1 %>%
    dplyr::select(-any_of("Age")) %>%
    merge(Part_Age %>% dplyr::select(response_uuid, Age),
          by = "response_uuid", all.x = TRUE) %>%
    distinct(response_uuid, Age) %>%
    group_by(Age) %>%
    summarise(n_usable = n(), .groups = "drop") %>%
    arrange(Age)
  
  sidebias_ids <- setdiff(ids_post_usable, ids_post_fixation)
  sidebias_by_age <- data.frame(response_uuid = sidebias_ids) %>%
    merge(Part_Age %>% dplyr::select(response_uuid, Age),
          by = "response_uuid", all.x = TRUE) %>%
    group_by(Age) %>%
    summarise(n_sidebias = n(), .groups = "drop") %>%
    merge(data.frame(Age = age_levels), all.y = TRUE) %>%
    mutate(n_sidebias = tidyr::replace_na(n_sidebias, 0)) %>%
    merge(usable_by_age, by = "Age") %>%
    mutate(pct_sidebias = sprintf("%.1f", n_sidebias / n_usable * 100)) %>%
    arrange(Age)
  
  final_by_age <- enrollment$metrics %>%
    dplyr::select(-any_of("Age")) %>%
    merge(Part_Age %>% dplyr::select(response_uuid, Age),
          by = "response_uuid", all.x = TRUE) %>%
    distinct(response_uuid, Age) %>%
    group_by(Age) %>%
    summarise(n_final = n(), .groups = "drop") %>%
    arrange(Age)
  
  list(
    task_label        = task_label,
    n_analytic_trials = n_analytic_trials,
    n_usable_1trial   = n_usable_1trial,
    n_sidebias        = n_sidebias,
    pct_sidebias      = pct_sidebias,
    final_n           = final_n,
    t_med = t_med, t_q1 = t_q1, t_q3 = t_q3,
    trials_by_age   = trials_by_age,
    usable_by_age   = usable_by_age,
    sidebias_by_age = sidebias_by_age,
    final_by_age    = final_by_age,
    age_levels      = age_levels
  )
}

make_table4_row <- function(enrollment, results, icatcher_df,
                            task_label, n_analytic_trials,
                            Part_Age, trial_col = "Trial",
                            final_df = NULL, final_id_col = "response_uuid",
                            final_trial_col = NULL,
                            trials_per_person_override = NULL) {
  
  age_levels  <- sort(unique(Part_Age$Age))
  age_counts  <- Part_Age %>%
    group_by(Age) %>%
    summarise(n_consented = n(), .groups = "drop") %>%
    arrange(Age)
  
  n_consented <- nrow(Part_Age)
  n_valid     <- length(unique(icatcher_df$response_uuid))
  pct_valid   <- round(n_valid / n_consented * 100, 1)
  
  # Valid iCatcher by age
  valid_by_age <- icatcher_df %>%
    merge(Part_Age %>% dplyr::select(response_uuid, Age),
          by = "response_uuid", all.x = TRUE) %>%
    distinct(response_uuid, Age) %>%
    group_by(Age) %>%
    summarise(n_valid = n(), .groups = "drop") %>%
    merge(age_counts, by = "Age") %>%
    mutate(pct_valid = sprintf("%.1f", n_valid / n_consented * 100)) %>%
    arrange(Age)
  
  row <- compute_table4_row(enrollment, task_label, n_analytic_trials,
                            Part_Age, trial_col)
  
  pct_usable <- round(row$n_usable_1trial / n_valid * 100, 1)
  
  # Usable by age (pct out of valid iCatcher per age)
  usable_by_age <- row$usable_by_age %>%
    merge(valid_by_age %>% dplyr::select(Age, n_valid), by = "Age") %>%
    mutate(pct_usable = sprintf("%.1f", n_usable / n_valid * 100)) %>%
    arrange(Age)
  
  sidebias_str <- if (row$n_sidebias == 0 &&
                      identical(sort(enrollment$metrics_1$response_uuid),
                                sort(enrollment$metrics$response_uuid))) {
    "---"
  } else {
    paste0(row$n_sidebias, " (", row$pct_sidebias, "%)")
  }
  
  # Final N and trials per person
  if (!is.null(final_df)) {
    
    t_col   <- if (!is.null(final_trial_col)) final_trial_col else trial_col
    final_n <- length(unique(final_df[[final_id_col]]))
    
    trials_per_person <- if (!is.null(trials_per_person_override)) {
      trials_per_person_override %>%
        dplyr::select(-any_of("Age")) %>%
        merge(Part_Age %>% dplyr::select(response_uuid, Age),
              by = "response_uuid", all.x = TRUE)
    } else {
      final_df %>%
        dplyr::select(-any_of("Age")) %>%
        merge(Part_Age %>% dplyr::select(response_uuid = !!final_id_col, Age),
              by.x = final_id_col, by.y = "response_uuid", all.x = TRUE) %>%
        group_by(.data[[final_id_col]], Age) %>%
        summarise(n_trials = n_distinct(.data[[t_col]]), .groups = "drop")
    }
    
    # Rounded to 1 dp
    t_med <- round(median(trials_per_person$n_trials), 1)
    t_q1  <- round(quantile(trials_per_person$n_trials, 0.25), 1)
    t_q3  <- round(quantile(trials_per_person$n_trials, 0.75), 1)
    
    trials_by_age <- trials_per_person %>%
      group_by(Age) %>%
      summarise(
        t_med = round(median(n_trials), 1),
        t_q1  = round(quantile(n_trials, 0.25), 1),
        t_q3  = round(quantile(n_trials, 0.75), 1),
        n     = n(),
        .groups = "drop"
      ) %>%
      arrange(Age)
    
    final_by_age <- final_df %>%
      dplyr::select(-any_of("Age")) %>%
      merge(Part_Age %>% dplyr::select(response_uuid = !!final_id_col, Age),
            by.x = final_id_col, by.y = "response_uuid", all.x = TRUE) %>%
      distinct(.data[[final_id_col]], Age) %>%
      group_by(Age) %>%
      summarise(n_final = n(), .groups = "drop") %>%
      arrange(Age)
    
  } else {
    final_n       <- row$final_n
    t_med         <- row$t_med
    t_q1          <- row$t_q1
    t_q3          <- row$t_q3
    trials_by_age <- row$trials_by_age
    final_by_age  <- row$final_by_age
  }
  
  # --- Build output row with split columns ---
  out <- data.frame(
    Task            = task_label,
    Analytic_trials = n_analytic_trials,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(age_levels)) {
    age <- age_levels[i]
    
    vba  <- valid_by_age        %>% filter(Age == age)
    uba  <- usable_by_age       %>% filter(Age == age)
    sba  <- row$sidebias_by_age %>% filter(Age == age)
    fba  <- final_by_age        %>% filter(Age == age)
    tba  <- trials_by_age       %>% filter(Age == age)
    
    out[[paste0("N_valid_",    age)]] <- paste0(vba$n_valid,  " (", vba$pct_valid,  "%)")
    out[[paste0("N_usable_",   age)]] <- paste0(uba$n_usable, " (", uba$pct_usable, "%)")
    out[[paste0("N_sidebias_", age)]] <- if (sba$n_sidebias == 0 && sidebias_str == "---") "---" else
      paste0(sba$n_sidebias, " (", sba$pct_sidebias, "%)")
    out[[paste0("N_final_",    age)]] <- fba$n_final
    out[[paste0("Med_trials_", age)]] <- paste0(
      sprintf("%.1f", tba$t_med), " (",
      sprintf("%.1f", tba$t_q1),  ", ",
      sprintf("%.1f", tba$t_q3),  ")"
    )
  }
  
  out
}

# ------------------------------------------------------------------------------
# make_summary_table()
#
# Compiles results from multiple run_condition_model_v2() calls into a single
# wide summary table with one row per task, grouped columns per age group,
# and an Overall column — matching the layout in the target figure.
#
# USAGE:
#   Each entry in `model_results` is a named list produced by a slightly
#   modified run_condition_model_v2() that returns the raw `tab` data frame
#   instead of (or in addition to) the kable. See note below.
#
#   results_list <- list(
#     list(label = "Geo-Social Attention", tab = tab_geo),
#     list(label = "VPC Social (face) trials", tab = tab_vpc_social),
#     ...
#   )
#   make_summary_table(results_list, age_levels = c(6, 12))
#
# NOTE on run_condition_model_v2():
#   Add `return_data = FALSE` parameter to that function, and when TRUE,
#   return `tab` as a data frame instead of rendering the kable. This lets
#   you collect `tab` for each task and pass it here.
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# make_summary_table()
#
# Compiles results from multiple run_condition_model() calls into a single
# wide summary table. Handles two tab formats:
#
#   Non-interaction tabs: Condition values are "Overall", "Age 12 - Age 6", "6", "12"
#   Interaction tabs:     Condition values are "Non-social | Overall",
#                         "Non-social | Age 6", "Social | Age 12", etc.
#                         These are exploded into one row per stimulus type.
#
# USAGE:
#   results_list <- list(
#     list(label = "Geo-Social Attention",      tab = tab_geo),        # non-interaction
#     list(label = "VPC",                        tab = tab_vpc,         # interaction
#                                                stim_levels = c("Non-social", "Social"),
#                                                stim_labels = c("VPC Non-social", "VPC Social")),
#     ...
#   )
#   make_summary_table(results_list, age_levels = c(6, 12))
#
# For interaction entries:
#   stim_levels  — the stimulus type strings as they appear in the Condition column
#                  (e.g. "Non-social", "Social")
#   stim_labels  — optional display labels for the Task column; defaults to
#                  paste(label, stim_level)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# make_summary_table()
#
# Compiles results from multiple run_condition_model() calls into a single
# wide summary table. Handles two tab formats:
#
#   Non-interaction tabs: Condition values are "Overall", "Age 12 - Age 6", "6", "12"
#   Interaction tabs:     Condition values are "Non-social | Overall",
#                         "Non-social | Age 6", "Social | Age 12", etc.
#                         These are exploded into one row per stimulus type.
#
# USAGE:
#   results_list <- list(
#     list(label = "Geo-Social Attention",      tab = tab_geo),        # non-interaction
#     list(label = "VPC",                        tab = tab_vpc,         # interaction
#                                                stim_levels = c("Non-social", "Social"),
#                                                stim_labels = c("VPC Non-social", "VPC Social")),
#     ...
#   )
#   make_summary_table(results_list, age_levels = c(6, 12))
#
# For interaction entries:
#   stim_levels  — the stimulus type strings as they appear in the Condition column
#                  (e.g. "Non-social", "Social")
#   stim_labels  — optional display labels for the Task column; defaults to
#                  paste(label, stim_level)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# make_summary_table()
#
# Compiles results from multiple run_condition_model() calls into a single
# wide summary table. Handles two tab formats:
#
#   Non-interaction tabs: Condition values are "Overall", "Age 12 - Age 6", "6", "12"
#   Interaction tabs:     Condition values are "Non-social | Overall",
#                         "Non-social | Age 6", "Social | Age 12", etc.
#                         These are exploded into one row per stimulus type.
#
# USAGE:
#   results_list <- list(
#     list(label = "Geo-Social Attention",      tab = tab_geo),        # non-interaction
#     list(label = "VPC",                        tab = tab_vpc,         # interaction
#                                                stim_levels = c("Non-social", "Social"),
#                                                stim_labels = c("VPC Non-social", "VPC Social")),
#     ...
#   )
#   make_summary_table(results_list, age_levels = c(6, 12))
#
# For interaction entries:
#   stim_levels  — the stimulus type strings as they appear in the Condition column
#                  (e.g. "Non-social", "Social")
#   stim_labels  — optional display labels for the Task column; defaults to
#                  paste(label, stim_level)
# ------------------------------------------------------------------------------

make_summary_table <- function(results_list,
                               age_levels,
                               caption = "Summary of mixed model estimates") {
  
  age_levels <- age_levels
  
  # ----------------------------------------------------------------------------
  # Helper: extract one wide row from a non-interaction tab
  # ----------------------------------------------------------------------------
  
  extract_non_interaction <- function(tab, label) {
    
    overall_row   <- tab[tab$Condition == "Overall",         , drop = FALSE]
    contrast_row  <- tab[tab$Condition == "Age 12 - Age 6",  , drop = FALSE]
    
    age_rows <- lapply(age_levels, function(a) {
      tab[tab$Condition == a, , drop = FALSE]
    })
    names(age_rows) <- age_levels
    
    out <- data.frame(Task = label, stringsAsFactors = FALSE)
    
    out[["Overall_Estimate"]]  <- if (nrow(overall_row))  overall_row[["Estimate (95% CI)"]]  else NA_character_
    out[["Overall_P"]]         <- if (nrow(overall_row))  overall_row[["P-value"]]             else NA_character_
    out[["Contrast_Estimate"]] <- if (nrow(contrast_row)) contrast_row[["Estimate (95% CI)"]] else NA_character_
    out[["Contrast_P"]]        <- if (nrow(contrast_row)) contrast_row[["P-value"]]            else NA_character_
    
    for (a in age_levels) {
      r <- age_rows[[a]]
      out[[paste0("Age_", a, "_Estimate")]] <- if (nrow(r)) r[["Estimate (95% CI)"]] else NA_character_
      out[[paste0("Age_", a, "_P")]]        <- if (nrow(r)) r[["P-value"]]            else NA_character_
    }
    
    out
  }
  
  # ----------------------------------------------------------------------------
  # Helper: extract one wide row per stimulus level from an interaction tab
  # Condition format: "<stim> | Overall", "<stim> | Age 6", "<stim> | Age 12"
  # ----------------------------------------------------------------------------
  
  extract_interaction <- function(tab, label, stim_levels, stim_labels = NULL) {
    
    if (is.null(stim_labels)) {
      stim_labels <- paste(label, stim_levels)
    }
    
    lapply(seq_along(stim_levels), function(i) {
      
      stim   <- stim_levels[i]
      rlabel <- stim_labels[i]
      
      overall_row  <- tab[tab$Condition == paste0(stim, " | Overall"),       , drop = FALSE]
      contrast_row <- tab[tab$Condition == paste0(stim, " | Age 12 - Age 6"), , drop = FALSE]
      
      age_rows <- lapply(age_levels, function(a) {
        tab[tab$Condition == paste0(stim, " | Age ", a), , drop = FALSE]
      })
      names(age_rows) <- age_levels
      
      out <- data.frame(Task = rlabel, stringsAsFactors = FALSE)
      
      out[["Overall_Estimate"]]  <- if (nrow(overall_row))  overall_row[["Estimate (95% CI)"]]  else NA_character_
      out[["Overall_P"]]         <- if (nrow(overall_row))  overall_row[["P-value"]]             else NA_character_
      out[["Contrast_Estimate"]] <- if (nrow(contrast_row)) contrast_row[["Estimate (95% CI)"]] else NA_character_
      out[["Contrast_P"]]        <- if (nrow(contrast_row)) contrast_row[["P-value"]]            else NA_character_
      
      for (a in age_levels) {
        r <- age_rows[[a]]
        out[[paste0("Age_", a, "_Estimate")]] <- if (nrow(r)) r[["Estimate (95% CI)"]] else NA_character_
        out[[paste0("Age_", a, "_P")]]        <- if (nrow(r)) r[["P-value"]]            else NA_character_
      }
      
      out
    })
  }
  
  # ----------------------------------------------------------------------------
  # Process each entry in results_list
  # ----------------------------------------------------------------------------
  
  rows <- lapply(results_list, function(res) {
    
    is_interaction <- !is.null(res$stim_levels)
    
    if (is_interaction) {
      extract_interaction(
        tab         = res$tab,
        label       = res$label,
        stim_levels = res$stim_levels,
        stim_labels = res$stim_labels   # NULL is fine; defaults to paste(label, stim)
      )
    } else {
      list(extract_non_interaction(res$tab, res$label))
    }
  })
  
  # Flatten the nested list (interaction entries produce multiple rows)
  wide <- do.call(rbind, do.call(c, rows))
  
  # ----------------------------------------------------------------------------
  # Reorder columns: Task | Overall | Age 6 | Age 12 | ...
  # Rename columns for display and build spanning header
  # ----------------------------------------------------------------------------
  
  age_col_pairs <- unlist(lapply(age_levels, function(a)
    c(paste0("Age_", a, "_Estimate"), paste0("Age_", a, "_P"))))
  
  wide <- wide[, c("Task",
                   "Overall_Estimate",  "Overall_P",
                   age_col_pairs,
                   "Contrast_Estimate", "Contrast_P"
  )]
  
  col_rename <- c("Task"              = "Task",
                  "Overall_Estimate"  = "Estimate (95% CI)",
                  "Overall_P"         = "P-value",
                  "Contrast_Estimate" = "Estimate (95% CI)",
                  "Contrast_P"        = "P-value")
  for (a in age_levels) {
    col_rename[paste0("Age_", a, "_Estimate")] <- "Estimate (95% CI)"
    col_rename[paste0("Age_", a, "_P")]        <- "P-value"
  }
  
  colnames(wide) <- col_rename
  
  header <- c(
    " "              = 1,
    "Overall"        = 2,
    setNames(rep(2, length(age_levels)), paste0("Age ", age_levels)),
    "Age 12 - Age 6" = 2
  )
  
  # ----------------------------------------------------------------------------
  # Render
  # ----------------------------------------------------------------------------
  
  wide %>%
    knitr::kable(
      format  = "html",
      caption = caption,
      align   = c("l", rep(c("c", "c"), length(age_levels) + 1))
    ) %>%
    kableExtra::kable_styling("striped", full_width = FALSE) %>%
    kableExtra::add_header_above(header)
}

make_summary_table_no_age <- function(results_list,
                                      caption = "Summary of mixed model estimates") {
  
  # ----------------------------------------------------------------------------
  # Helper: single-row entry (no stim levels)
  # ----------------------------------------------------------------------------
  
  extract_single <- function(tab, label) {
    row <- tab[1, , drop = FALSE]
    data.frame(
      Task             = label,
      Overall_Estimate = row[["Estimate (95% CI)"]],
      Overall_P        = row[["P-value"]],
      stringsAsFactors = FALSE
    )
  }
  
  # ----------------------------------------------------------------------------
  # Helper: one row per stimulus level
  # ----------------------------------------------------------------------------
  
  extract_stim <- function(tab, label, stim_levels, stim_labels = NULL) {
    if (is.null(stim_labels)) stim_labels <- paste(label, stim_levels)
    
    lapply(seq_along(stim_levels), function(i) {
      stim <- stim_levels[i]
      row  <- tab[tab$Condition == stim, , drop = FALSE]
      data.frame(
        Task             = stim_labels[i],
        Overall_Estimate = if (nrow(row)) row[["Estimate (95% CI)"]] else NA_character_,
        Overall_P        = if (nrow(row)) row[["P-value"]]           else NA_character_,
        stringsAsFactors = FALSE
      )
    })
  }
  
  # ----------------------------------------------------------------------------
  # Process each entry
  # ----------------------------------------------------------------------------
  
  rows <- lapply(results_list, function(res) {
    if (!is.null(res$stim_levels)) {
      extract_stim(
        tab         = res$tab,
        label       = res$label,
        stim_levels = res$stim_levels,
        stim_labels = res$stim_labels
      )
    } else {
      list(extract_single(res$tab, res$label))
    }
  })
  
  wide <- do.call(rbind, do.call(c, rows))
  
  # ----------------------------------------------------------------------------
  # Rename and render
  # ----------------------------------------------------------------------------
  
  colnames(wide) <- c("Task", "Estimate (95% CI)", "P-value")
  
  wide %>%
    knitr::kable(
      format  = "html",
      caption = caption,
      align   = c("l", "c", "c")
    ) %>%
    kableExtra::kable_styling("striped", full_width = FALSE) %>%
    kableExtra::add_header_above(c(" " = 1, "Overall" = 2))
}
