# ------------------------------------------------------------
# compute_gaze_metrics()
#
# Creates QC-filtered analytic gaze summaries from
# iCatcher processing outputs.
#
# Supports:
#   - dual-screen tasks
#   - single-screen tasks
#   - VPC/HelperHinderer trial collapsing
#
# Main steps:
#   1. Remove low usable-frame trials
#   2. Optionally remove fixation participants
#   3. Compute trial-level gaze metrics
#   4. Optionally collapse VPC/HelperHinderer trials
#   5. Optionally collapse VPC/HelperHinderer trials
#   6. Compute participant-level averages

# Inputs:
#
#   results:
#     output from process_icatcher() or
#     process_icatcher_single()
#
#   Part_Age:
#     participant age lookup table
#
#   task_type:
#     "dual" or "single"
#
#   usable_threshold:
#     minimum usable-frame proportion
#
#   remove_fixations:
#     whether to exclude fixation participants
#
#   remove_extreme_target:
#     whether to exclude trials where target_prop is 0 or 1
#     (dual-screen only)

#   collapse_vpc:
#     whether to average across VPC/ (LRISE) helperhinderer trial pairs
#
#   type_task:
#     VPC/ (LRISE) helperhinderer task lookup table (Social vs Non-social)
#
#   task_icatcher:
#     metadata table containing Version
#
#   trial_filter:
#     optional character/numeric vector of Trial values to keep;
#     applied before all other steps. NULL = keep all trials.

#
# Returns:
#
#   list containing:
#     - metrics_1   : post-usable-threshold, pre-fixation removal
#     - metrics     : final trial-level analytic sample
#     - metrics_avg : participant-level averages
# ------------------------------------------------------------

compute_gaze_metrics <- function(results,
                                 Part_Age,
                                 task_type,
                                 usable_threshold,
                                 remove_fixations,
                                 remove_extreme_target,
                                 collapse_vpc_hh,
                                 type_task      = NULL,
                                 task_icatcher  = NULL,
                                 trial_filter   = NULL) {   # NULL = keep all trials
  
  task_type <- match.arg(task_type, choices = c("dual", "single"))
  
  # Filter to specific trials if requested
  if (!is.null(trial_filter)) {
    if (!is.null(results$side_bytrial)) 
      results$side_bytrial <- results$side_bytrial %>% filter(Trial %in% trial_filter)
    if (!is.null(results$trial_counts))
      results$trial_counts <- results$trial_counts %>% filter(Trial %in% trial_filter)
    if (!is.null(results$target_bytrial))
      results$target_bytrial <- results$target_bytrial %>% filter(Trial %in% trial_filter)
    if (!is.null(results$quality_flags))
      results$quality_flags <- results$quality_flags %>% filter(Trial %in% trial_filter)
  }
  
  # ----------------------------------------------------------
  # Merge trial counts with quality flags and compute
  # iCatcher coverage and usable frame proportion (single screen)
  # ----------------------------------------------------------
  
  if (task_type == "dual") {
    
    df_metrics_1 <- merge(
      results$trial_counts,
      results$target_bytrial,
      all = TRUE
    )
    
  } else {
    
    df_metrics_1 <- results$trial_counts %>%
      merge(results$quality_flags %>% dplyr::select(-Age),
            by  = c("response_uuid", "Trial"),
            all = TRUE)%>%
      dplyr::mutate(
        has_icatcher  = ifelse(noicatcherdata != frames | is.na(noicatcherdata), T, F) 
      )%>%
      mutate(noicatcherdata = tidyr::replace_na(noicatcherdata, 0),
             n_icatcher     = frames - noicatcherdata) %>%
      dplyr::select(response_uuid, Trial, count_85_non_na, n_icatcher, frames)
    
  }
  
  # ----------------------------------------------------------
  # Compute usable frame proportion
  # ----------------------------------------------------------
  
  df_metrics_1 <- df_metrics_1 %>%
    mutate(
      usable_prop =
        round(
          count_85_non_na / frames,
          3
        )
    ) %>%
    merge(Part_Age, all.x = TRUE) %>%
    filter(
      usable_prop >= usable_threshold
    )
  
  # ----------------------------------------------------------
  # Remove fixation participants
  # Dual-screen tasks only
  # ----------------------------------------------------------
  
  if (task_type == "dual" && remove_fixations) {
    
    fixation_summary <- results$side_bytrial %>%
      merge(
        df_metrics_1 %>%
          dplyr::select(
            response_uuid,
            Trial
          )
      ) %>%
      group_by(response_uuid) %>%
      mutate(
        left = sum(left),
        right = sum(right)
      ) %>%
      dplyr::select(-Trial) %>%
      distinct() %>%
      mutate(
        total_looks = left + right,
        left_prop = left / total_looks,
        right_prop = right / total_looks,
        fixation_flag =
          ifelse(
            left_prop >= .80 |
              right_prop >= .80,
            1,
            0
          )
      ) %>%
      filter(fixation_flag == 0) %>%
      ungroup()
    
    df_metrics <- df_metrics_1 %>%
      filter(
        response_uuid %in%
          fixation_summary$response_uuid
      )
    
  } else {
    
    df_metrics <- df_metrics_1
    
  }
  
  # ----------------------------------------------------------
  # Compute task-specific gaze metric:
  #   dual   — target_prop = target / (target + non_target)
  #   single — onscreen_prop = high-confidence frames / iCatcher frames
  # ----------------------------------------------------------
  
  if (task_type == "dual") {
    
    df_metrics <- df_metrics %>%
      mutate(
        target_prop =
          round(
            target / (non_target + target),
            3
          )
      )
    
    if (remove_extreme_target) {
      
      df_metrics <- df_metrics %>%
        filter(
          !is.nan(target_prop),
          !target_prop %in% c(0, 1)
        )
      
    } else {
      
      df_metrics <- df_metrics %>%
        filter(
          !is.nan(target_prop)
        )
      
    }
    
  } else {
    
    df_metrics <- df_metrics %>%
      mutate(
        onscreen_prop =
          round(
            (count_85_non_na) /
              (n_icatcher),
            3
          )
      )
    
  }
  
  # ----------------------------------------------------------
  # Collapse VPC/HelperHinderer trial pairs
  # ----------------------------------------------------------
  
  if (collapse_vpc_hh) {
    
    if (is.null(type_task) || is.null(task_icatcher)) {
      stop("collapse_vpc_hh = TRUE requires type_task and task_icatcher to be provided.")
    }
    
    # Build the trial-pair lookup: Trial (character) -> Trial Total + Stimulus type
    trial_lookup <- type_task %>%
      mutate(
        `Trial Total` = as.numeric(gsub("[^0-9]", "", `Trial Number`)),
        Trial = as.character(`Trial Number`)       # ensure type matches df_metrics$Trial
      ) %>%
      dplyr::select(Trial, `Trial Total`, `Stimulus type`)
    
    # Attach Version from task_icatcher
    version_lookup <- task_icatcher %>%
      dplyr::select(response_uuid, Version) %>%
      distinct()
    
    df_metrics <- df_metrics %>%
      merge(version_lookup, by = "response_uuid") %>%
      merge(trial_lookup,   by = "Trial") %>%
      group_by(response_uuid, `Trial Total`) %>%
      mutate(
        target_prop = mean(target_prop, na.rm = TRUE)  # average across sub-trials in pair
      ) %>%
      arrange(response_uuid, `Trial Total`) %>%
      slice(1) %>%
      ungroup() %>%
      dplyr::select(-Trial)                            # drop sub-trial; Trial Total is the unit now
    
    # Verify the merge worked — warn if any rows were lost unexpectedly
    if (nrow(df_metrics) == 0) {
      warning("collapse_vpc_hh merge produced 0 rows. Check that Trial values in df_metrics match Trial Number in type_task.")
    }
  }
  
  # ----------------------------------------------------------
  # Compute participant-level averages
  # ----------------------------------------------------------
  
  metric_col <- ifelse(
    task_type == "dual",
    "target_prop",
    "onscreen_prop"
  )
  
  metric_avg_col <- ifelse(
    task_type == "dual",
    "target_prop_avg",
    "onscreen_prop_avg"
  )
  
  df_metrics_avg <- df_metrics %>%
    group_by(response_uuid) %>%
    summarise(
      metric_avg =
        mean(.data[[metric_col]]),
      .groups = "drop"
    ) %>%
    rename(
      !!metric_avg_col := metric_avg
    ) %>%
    merge(
      Part_Age,
      all.x = TRUE
    )
  
  # ----------------------------------------------------------
  # Return outputs
  # ----------------------------------------------------------
  
  list(
    metrics_1   = df_metrics_1,
    metrics     = df_metrics,
    metrics_avg = df_metrics_avg
  )
}
