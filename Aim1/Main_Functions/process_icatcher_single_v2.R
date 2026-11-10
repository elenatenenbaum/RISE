# ------------------------------------------------------------
# process_icatcher_single()
#
# Variant of process_icatcher() for single-screen tasks where
# there is no target position — only left/right side looking
# is recorded, not target vs non-target.
#
# Main steps:
#   1. Load and clean raw iCatcher output files
#   2. Remove leading frames specified in metadata
#   3. Merge frame labels and participant age information
#   4. Compute quality-control summaries:
#        - look away frames
#        - no-face frames
#        - low-confidence frames
#        - missing iCatcher frames
#   5. Compute trial-level valid frame counts and proportions
#   6. Classify high-confidence looks as left or right only
#        (no target/non-target — no position column)
#   7. Return frame-level and trial-level summary outputs
#
# Inputs:
#   icatcher_data  : metadata table containing file paths,
#                    response_uuid, frame trimming info,
#                    and optional Version column
#
#   labeled_frames : frame-by-frame task labels containing
#                    Trial information (no position column)
#
#   Part_Age       : participant age lookup table
#
# Returns:
#   A named list containing:
#     - total_rows
#     - frame_counts
#     - quality_flags
#     - trial_counts
#     - side_bytrial
#     - frame_level
#   (target_bytrial omitted — no position in single-screen task)
# ------------------------------------------------------------

process_icatcher_single <- function(
    icatcher_data,
    labeled_frames,
    Part_Age,
    use_version = FALSE
) {
  
  icatcher_files <- icatcher_data$icatcher_files
  
  # ----------------------------------------------------------
  # Initialize output containers
  # ----------------------------------------------------------
  
  total_rows_list      <- list()
  frame_counts_list    <- list()
  quality_flags_list   <- list()
  trial_counts_list    <- list()
  side_bytrial_list    <- list()
  frame_level_list     <- list()
  
  # ----------------------------------------------------------
  # Process each iCatcher file
  # ----------------------------------------------------------
  
  for (file in icatcher_files) {
    
    participant_id <- icatcher_data %>%
      filter(icatcher_files == file) %>%
      pull(response_uuid) %>%
      unique()
    
    stopifnot(length(participant_id) == 1)
    
    # --------------------------------------------------------
    # Load and clean raw iCatcher output
    # --------------------------------------------------------
    
    merge_cols <- c("response_uuid", "frame_to_remove")
    
    if (use_version) {
      merge_cols <- c(merge_cols, "Version")
    }
    
    raw_icatcher <- read.csv(file, header = FALSE) %>%
      mutate_all(str_trim) %>%
      rename(frame = V1) %>%
      mutate(
        frame         = as.numeric(frame),
        response_uuid = participant_id
      ) %>%
      merge(
        icatcher_data %>%
          dplyr::select(all_of(merge_cols))
      )
    
    # --------------------------------------------------------
    # Remove leading frames
    # --------------------------------------------------------
    
    n_remove <- raw_icatcher$frame_to_remove[1]
    
    if (n_remove > 0) {
      trimmed_icatcher <- raw_icatcher %>%
        slice(-(1:n_remove))
    } else {
      trimmed_icatcher <- raw_icatcher
    }
    
    trimmed_icatcher <- trimmed_icatcher %>%
      mutate(frame = row_number() - 1)
    
    total_rows_list[[participant_id]] <- trimmed_icatcher %>%
      summarise(n = n())
    
    # --------------------------------------------------------
    # Select labeled frames
    # For single screen: version selects Trial column only,
    # no position column involved
    # --------------------------------------------------------
    
    if (use_version) {
      
      version_col <- unique(trimmed_icatcher$Version)
      
      stopifnot(length(version_col) == 1)
      
      version_col <- version_col[1]
      
      if (!version_col %in% names(labeled_frames)) {
        stop(
          paste(
            "Version column",
            version_col,
            "not found in labeled_frames"
          )
        )
      }
      
      # Single screen: version determines Trial column only,
      # no position selected
      labeled_subset <- labeled_frames %>%
        filter(!is.na(get(version_col)))%>%
        dplyr::select(
          frame,
          Trial
        ) 
      
    } else {
      
      labeled_subset <- labeled_frames
      
    }
    
    # --------------------------------------------------------
    # Merge iCatcher output with trial labels and age data
    # No position column for single-screen task
    # --------------------------------------------------------
    
    merged_data <- merge(
      trimmed_icatcher,
      labeled_subset,
      all.y = TRUE,
      by    = "frame"
    ) %>%
      # Fill response_uuid first so Part_Age merge can match
      mutate(response_uuid = participant_id) %>%  # restore if all.y merge drops it
      merge(Part_Age, all.x = TRUE) %>%
      tidyr::fill(Age, .direction = "downup") %>%
      arrange(frame)
    
    # --------------------------------------------------------
    # Compute quality-control flags
    # --------------------------------------------------------
    
    lookaway_summary <- merged_data %>%
      filter(V2 == "away") %>%
      group_by(Age, Trial) %>%
      summarise(away = n(), .groups = "drop")
    
    noface_summary <- merged_data %>%
      filter(V2 %in% c("nobabyface", "noface", "none")) %>%
      group_by(Age, Trial) %>%
      summarise(noface = n(), .groups = "drop")
    
    lowconf_summary <- merged_data %>%
      filter(V2 %in% c("left", "right"), V3 < 0.85) %>%
      group_by(Age, Trial) %>%
      summarise(lowconf = n(), .groups = "drop")
    
    missingdata_summary <- merged_data %>%
      filter(is.na(V2)) %>%
      group_by(Age, Trial) %>%
      summarise(noicatcherdata = n(), .groups = "drop")
    
    quality_flags <- lookaway_summary %>%
      full_join(noface_summary,      by = c("Age", "Trial")) %>%
      full_join(lowconf_summary,     by = c("Age", "Trial")) %>%
      full_join(missingdata_summary, by = c("Age", "Trial"))
    
    quality_flags_list[[participant_id]] <- quality_flags
    
    # --------------------------------------------------------
    # Count raw frames per trial (all frames exist in task video, pre-confidence filtering)
    # --------------------------------------------------------
    
    frame_counts <- merged_data %>%
      group_by(Trial) %>%
      summarise(n = n(), .groups = "drop")
    
    frame_counts_list[[participant_id]] <- frame_counts
    
    # --------------------------------------------------------
    # Restrict to valid left/right looks
    # --------------------------------------------------------
    
    merged_data <- merged_data %>%
      mutate(
        V3 = as.numeric(V3),
        V2 = if_else(
          V2 %in% c("left", "right"),
          V2,
          NA_character_
        ),
        V3 = if_else(is.na(V2), NA_real_, V3)
      )
    
    # --------------------------------------------------------
    # Compute valid-frame summaries
    # Single screen:
    #   - frames = all rows (denominator)
    #   - count_85_non_na = left or right at >= 85% (numerator)
    #   - prop_looking_85 = count_85_non_na / frames
    # --------------------------------------------------------
    
    trial_counts <- merged_data %>%
      group_by(Trial) %>%
      summarise(
        non_na_count    = sum(!is.na(V2)),               # any left/right regardless of conf
        frames          = n(),                            # all rows = denominator
        count_85_non_na = sum(!is.na(V2) & V3 >= 0.85), # left or right at >= 85% (numerator)
        .groups = "drop"
      )
    
    trial_counts_list[[participant_id]] <- trial_counts
    
    # --------------------------------------------------------
    # Restrict to high-confidence looks
    # No target/non-target classification for single screen
    # --------------------------------------------------------
    
    highconf_data <- merged_data %>%
      filter(V3 >= 0.85)
    
    # --------------------------------------------------------
    # Compute side summary (left vs right at >= 85%)
    # target_bytrial removed — no position in single screen
    # --------------------------------------------------------
    
    if (nrow(highconf_data) == 0) {
      
      side_bytrial <- data.frame(
        Var1  = NA,
        left  = NA,
        right = NA
      )
      
    } else {
      
      side_bytrial <- as.data.frame(
        table(highconf_data$Trial, highconf_data$V2)
      ) %>%
        tidyr::spread(Var2, Freq) %>%
        { if (!"left"  %in% names(.)) mutate(., left  = 0L) else . } %>%
        { if (!"right" %in% names(.)) mutate(., right = 0L) else . }
      
    }
    
    side_bytrial_list[[participant_id]] <- side_bytrial
    
    # --------------------------------------------------------
    # Create frame-level output
    # look: "left"/"right" at >= 85%, NA otherwise
    # position removed — no position in single screen
    # --------------------------------------------------------
    
    frame_level <- merged_data %>%
      mutate(
        look = if_else(
          V3 >= 0.85 & !is.na(V2),
          V2,            # "left" or "right"
          NA_character_  # away, noface, low conf, missing
        )
      ) %>%
      dplyr::select(
        response_uuid,
        Trial,
        frame,
        V2,
        V3,
        look
      ) %>%
      group_by(Trial) %>%
      mutate(
        frame_within_trial = row_number()
      ) %>%
      ungroup()
    
    frame_level_list[[participant_id]] <- frame_level
  }
  
  # ----------------------------------------------------------
  # Combine participant-level outputs
  # ----------------------------------------------------------
  
  df_total_rows <- bind_rows(
    total_rows_list,
    .id = "response_uuid"
  )
  
  df_frame_counts <- bind_rows(
    frame_counts_list,
    .id = "response_uuid"
  )
  
  df_quality_flags <- bind_rows(
    quality_flags_list,
    .id = "response_uuid"
  )
  
  df_trial_counts <- bind_rows(
    trial_counts_list,
    .id = "response_uuid"
  )
  
  df_side_bytrial <- bind_rows(
    side_bytrial_list,
    .id = "response_uuid"
  ) %>%
    rename(Trial = Var1)
  
  df_frame_level <- bind_rows(frame_level_list)
  
  # ----------------------------------------------------------
  # Return outputs
  # target_bytrial removed — no position in single screen task
  # ----------------------------------------------------------
  
  list(
    total_rows     = df_total_rows,
    frame_counts   = df_frame_counts,
    quality_flags  = df_quality_flags,
    trial_counts   = df_trial_counts,
    side_bytrial   = df_side_bytrial,
    frame_level    = df_frame_level
  )
}