# ============================================================
# Find long videos - Aim 2
# ============================================================

library(readxl)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)

# ---- File paths (edit these two lines only) -----------------------------
main_dir   <-  here::here("Aim2") 
output_dir <- here::here("Aim2", "Data", "Long_Video") 
folder_path <- file.path(main_dir,"Data", "Annotation")
base_dir    <- file.path(main_dir,"Data", "Frame_Data")

# Expected full length (in frames) for each video type --------------------
# used to compute what proportion of the video was actually recorded.
folder_path_video <- file.path(main_dir, "Data", "Video")

files_video <- list.files(folder_path_video, pattern = "\\.(mp4|mov|avi|mkv)$",
                          full.names = TRUE, ignore.case = TRUE)

get_duration <- function(f) {
  out <- system2("ffprobe",
                 args = c("-v", "error",
                          "-show_entries", "format=duration",
                          "-of", "default=noprint_wrappers=1:nokey=1",
                          shQuote(f)),
                 stdout = TRUE)
  as.numeric(out)
}

durations <- sapply(files_video, get_duration)
frames <- round(durations * 30.3)
results <- data.frame(file = basename(files_video), seconds = durations, frames = frames)

video_lengths <- results %>%
  mutate(task = case_when(
    grepl("^Geo",   file) ~ "Geo",
    grepl("^AV", file)          ~ "AV",
    grepl("^wordrecog",   file) ~ "word",
    grepl("^VPCV",        file) ~ "VPC",
    grepl("^Prediction",  file) ~ "Pred",
    grepl("^Numeracy",    file) ~ "Num",
    grepl("^condition",   file) ~ "condi"
  )) %>%
  group_by(task) %>%
  summarise(frames = round(median(frames))) %>%
  { setNames(.$frames, .$task) }


# Expected full length (in frames) for each video type --------------------
# used to compute what proportion of the video was actually recorded.
# video_lengths <- c(
#   Geo   = 1712,
#   AV    = 3062,
#   word  = 4355,
#   VPC   = 8949,
#   Pred  = 1952,
#   Num   = 3213,
#   condi = 10661
# )

# Vectorized helper: looks up the matching video_lengths entry for each
# element of `value`, and returns count / expected_length.
# (Replaces the two duplicated case_when blocks that used to live below.)

read_session_folder <- function(subfolder, fix_event_number = FALSE) {
  files <- list.files(file.path(base_dir, subfolder), full.names = TRUE)
  df_list <- lapply(files, read_csv)
  
  if (fix_event_number) {
    # event_number is read as character in some files here; coerce before binding
    df_list <- lapply(df_list, function(d) {
      as.data.frame(d) %>% dplyr::mutate(event_number = as.numeric(event_number))
    })
  }
  
  bind_rows(df_list)
}

calc_prop <- function(value, count) {
  denom <- rep(NA_real_, length(value))
  for (pattern in names(video_lengths)) {
    denom[grepl(pattern, value)] <- video_lengths[[pattern]]
  }
  count / denom
}

# ---- Import annotated (iCatcher) data ------------------------------------
# takes only test trials
icatcher_files <- list.files(
  path = folder_path,
  recursive = TRUE,
  #pattern = "*test-trial*|*test-post*",
  full.names = TRUE
)
icatcher_file_names <- basename(icatcher_files)

# frame_id / response_uuid are the first two "_"-separated chunks of the file name
icatcher_output <- data.frame(
  icatcher_files      = icatcher_files,
  icatcher_file_names = icatcher_file_names,
  frame_id            = word(icatcher_file_names, 1, sep = "_"),
  response_uuid       = word(icatcher_file_names, 2, sep = "_")
)

# TODO: remove id in babysibs, in retest, and extra IDs in consent()

# ---- Import raw frame data -----------------------------------------------
raw_data_1 <- read_session_folder("RISE-Study_framedata_per_session")
raw_data_2 <- read_session_folder("RISE-Study-for-SPARK_framedata_per_session")
raw_data_3 <- read_session_folder("RISE-Study-Surveys_framedata_per_session", fix_event_number = TRUE)
raw_data_4 <- read_session_folder("RISE-Study-Follow-Up-Videos_framedata_per_session")

raw_data <- bind_rows(raw_data_1, raw_data_2, raw_data_3, raw_data_4) %>%
  distinct()

Stim_Order <- raw_data[raw_data$key == "videoShown", ]
# reduce to rows that contained "test"
Stim_Order <- Stim_Order[grep("test", Stim_Order$frame_id), ]

# merge participant order with icatcher output
icatcher <- merge(
  x = Stim_Order, y = icatcher_output,
  by = c("response_uuid", "frame_id"), all = FALSE
)

# ---- Count frames (rows) per video / icatcher file -----------------------
frame_ct <- list()
for (file in icatcher$icatcher_files) {
  id        <- as.character(icatcher$response_uuid[icatcher$icatcher_files == file])
  df        <- read.csv(file, header = FALSE)
  file_name <- icatcher$icatcher_file_names[icatcher$icatcher_files == file]
  
  frame_ct[[id]][[file_name]] <- nrow(df)
}

fr_ct <- bind_rows(frame_ct, .id = "response_uuid") %>%
  tidyr::pivot_longer(
    cols = -response_uuid,
    names_to = "icatcher_file_names",
    values_to = "row_count"
  ) %>%
  filter(!is.na(row_count)) %>%
  merge(icatcher %>% dplyr::select(1, value, icatcher_file_names))

df_ct <- fr_ct %>%
  mutate(prop = calc_prop(value, row_count)) %>%
  filter(prop >= 1.04) %>%
  arrange(prop) %>%
  mutate(prop = round(prop, 2))

hist(df_ct$prop, breaks = 16)

#View(fr_ct %>% group_by(value) %>% summarise(median = median(row_count)))
# A TOTAL OF 39 VIDEOS had proportion >= 1.05 (so 5% or more longer)
# If remove videos 10% or more longer: 34 videos
dim(df_ct %>% filter(prop >= 1.1))
dim(df_ct %>% filter(prop >= 1.2))
dim(df_ct %>% filter(prop >= 1.3))

write.csv(
  df_ct,
  file.path(output_dir, "aim2_total_frame_count_longvideo_beforetrim.csv"),
  row.names = FALSE
)

# ---- Pull timestamps from frame data --------------------------------------
video_time <- raw_data[grepl("test", raw_data$frame_id, ignore.case = TRUE), ]

options(digits.secs = 3)

video_time <- video_time %>%
  filter(key %in% c("eventType", "timestamp")) %>%
  tidyr::pivot_wider(names_from = key, values_from = value) %>%
  dplyr::select(response_uuid, child_hashed_id, frame_id, event_number, eventType, timestamp) %>%
  tidyr::pivot_wider(names_from = eventType, values_from = timestamp) %>%
  mutate(across(
    starts_with("exp-lookit"),
    ~ as.POSIXct(.x, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
  ))

# timestamps of startRecording, videoStarted, videoPaused
check <- video_time %>%
  tidyr::pivot_longer(
    cols = starts_with("exp-lookit-video"),
    names_to = "event_type",
    values_to = "timestamp"
  ) %>%
  filter(!is.na(timestamp)) %>%
  dplyr::select(1:4, 13:14)

clean_df <- video_time %>%
  dplyr::select(
    1:3,
    `exp-lookit-video:startRecording`,
    `exp-lookit-video:videoStarted`,
    `exp-lookit-video:videoPaused`,
    `exp-lookit-video:trialCompleted`
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("exp-lookit-video"),
    names_to = "event_type",
    values_to = "timestamp"
  ) %>%
  filter(!is.na(timestamp)) %>%
  group_by(response_uuid, child_hashed_id, frame_id, event_type) %>%
  arrange(timestamp, .by_group = TRUE) %>%
  mutate(event_idx = row_number()) %>%
  ungroup() %>%
  mutate(event_type = if_else(
    event_idx == 1,
    event_type,
    paste0(event_type, "_", event_idx)
  )) %>%
  dplyr::select(response_uuid, child_hashed_id, frame_id, event_type, timestamp) %>%
  tidyr::pivot_wider(names_from = event_type, values_from = timestamp)

clean_df <- clean_df %>%
  dplyr::select(
    1:3,
    `exp-lookit-video:startRecording`,
    `exp-lookit-video:videoStarted`,
    `exp-lookit-video:videoPaused`,
    `exp-lookit-video:videoStarted_2`,
    `exp-lookit-video:videoPaused_2`,
    `exp-lookit-video:trialCompleted`
  )

# time (secs) between recording start and video start, and pause/restart
df_diff <- clean_df %>%
  mutate(record_1ststart = as.numeric(difftime(
    `exp-lookit-video:videoStarted`, `exp-lookit-video:startRecording`, units = "secs"
  )))
# mutate(pause1_start2 = as.numeric(difftime(
#   `exp-lookit-video:videoPaused`, `exp-lookit-video:trialCompleted`, units = "secs"
# )))

df_ct <- df_ct %>%
  mutate(frame_id = sub("_.*$", "", icatcher_file_names))

# ---- Remove extra frames and recalculate proportions ----------------------
# converts seconds to number of frames/rows and recalculates the new proportions
new_df <- df_ct %>%
  merge(df_diff) %>%
  arrange(prop) %>%
  mutate(frame_to_remove = round((record_1ststart * 1000) / 33.3)) %>%
  mutate(frame_to_remove = ifelse(is.na(frame_to_remove), 0, frame_to_remove)) %>%
  mutate(new_rcount = row_count - frame_to_remove) %>%
  mutate(new_prop = calc_prop(value, new_rcount)) %>%
  filter(new_prop >= 1.04) %>%
  arrange(new_prop) %>%
  mutate(new_prop = round(new_prop, 2))

new_df <- new_df %>%
  dplyr::select(1:5, new_rcount, 14, 15, 17) %>%
  mutate(record_1ststart = round(record_1ststart, 1))

write.csv(
  new_df %>% arrange(response_uuid),
  file.path(output_dir, "aim2_total_frame_count_longvideo.csv"),
  row.names = FALSE
)