library(dplyr)
library(tidyr)
library(lubridate)
library(circular)
library(purrr)
library(ggplot2)

source("helpers.R")

#read files
ios_abn     <- read.csv("va_new_ios_logs.csv")
android_abn <- read.csv("va_new_android_logs.csv")

#filter for VAHR (easier to separate VAHR and VAPT into two files)
ios_vahr <- ios_abn %>%
  filter(grepl("^VAHR", association_code))

android_vahr <- android_abn %>%
  filter(grepl("^VAHR", association_code))

ios_clean <- ios_vahr %>%
  transmute(
    user_id    = users_idx,
    start_time = ymd_hms(start),
    end_time   = ymd_hms(end),
    duration   = as.numeric(duration),
    local_date = as.Date(local_date),
    raw_cat    = category,
    platform   = "ios"
  )

android_clean <- android_vahr %>%
  transmute(
    user_id    = users_idx,
    start_time = ymd_hms(start),
    end_time   = ymd_hms(end),
    duration   = as.numeric(duration),
    local_date = as.Date(local_date),
    raw_cat    = genre,
    platform   = "android"
  )

sessions <- bind_rows(ios_clean, android_clean) %>%
  mutate(
    start_hour = hour(start_time) + minute(start_time)/60,
    raw_cat    = ifelse(raw_cat == "" | is.na(raw_cat), "Unknown", raw_cat)
  ) %>%
  left_join(category_map, by = "raw_cat") %>%
  mutate(
    harmonized_category = if_else(is.na(harmonized_category),
                                  raw_cat,
                                  harmonized_category)
  )

#daily metrics
daily <- sessions %>%
  group_by(user_id, local_date) %>%
  summarise(
    total_dur      = sum(duration),
    n_sessions     = n(),
    mean_sess_dur  = mean(duration),
    mean_hour      = mean(start_hour),
    entropy        = shannon_entropy(harmonized_category),
    .groups = "drop"
  )

daily_center <- sessions %>%
  group_by(user_id, local_date) %>%
  summarise(
    center_time = sum(start_hour * duration, na.rm = TRUE) /
      sum(duration, na.rm = TRUE),
    .groups = "drop"
  )

daily <- daily %>%
  left_join(daily_center, by = c("user_id", "local_date")) %>%
  group_by(user_id) %>%
  arrange(local_date) %>%
  mutate(day_index = row_number()) %>%
  ungroup()

#overall variability calculations
eligible_users <- daily %>%
  count(user_id) %>%
  filter(n >= 21) %>% #can change n
  pull(user_id)

daily_eligible <- daily %>%
  filter(user_id %in% eligible_users)

variability_overall_vahr <- daily_eligible %>%
  group_by(user_id) %>%
  summarise(
    VEV_raw = 0.5 * (sd(total_dur)  / mean(total_dur)) +
      0.5 * (sd(n_sessions) / mean(n_sessions)),
    SSV_raw = sd(mean_sess_dur) / mean(mean_sess_dur),
    .groups = "drop"
  )

TCV_all <- daily_eligible %>%
  group_by(user_id) %>%
  summarise(
    TCV_raw = {
      theta <- circular(center_time * 2*pi/24, type = "angles", units = "radians")
      sd.circular(theta)
    },
    .groups = "drop"
  )

daily_cat <- sessions %>%
  filter(user_id %in% eligible_users) %>%
  group_by(user_id, local_date, harmonized_category) %>%
  summarise(total = sum(duration), .groups = "drop") %>%
  group_by(user_id, local_date) %>%
  mutate(p = total / sum(total)) %>%
  select(user_id, local_date, harmonized_category, p) %>%
  pivot_wider(names_from = harmonized_category,
              values_from = p,
              values_fill = 0) %>%
  arrange(user_id, local_date)

CCV_all <- daily_cat %>%
  group_by(user_id) %>%
  summarise(
    CCV_raw = {
      mats <- as.matrix(cur_data()[, -c(1,2)])
      if (nrow(mats) < 2) NA_real_ else {
        js_vals <- map_dbl(2:nrow(mats), ~ jsd(mats[.x-1, ], mats[.x, ]))
        mean(js_vals)
      }
    },
    .groups = "drop"
  )

variability_overall_vahr <- variability_overall_vahr %>%
  left_join(TCV_all, by = "user_id") %>%
  left_join(CCV_all, by = "user_id") %>%
  mutate(cohort = "VAHR")

#write data for external analysis
write.csv(variability_overall_vahr,
          "vahr_overall_variability_raw.csv",
          row.names = FALSE)

#sliding window analysis
window_size <- 21
all_window_rows <- list()

#loop through every eligible user and calculate BVS
for (u in eligible_users) {
  message("Processing user: ", u)
  windows <- get_user_windows(daily, u, window_size = window_size)
  if (length(windows) == 0) next
  
  for (w_name in names(windows)) {
    message(u, w_name)
    w_df <- windows[[w_name]]
    scores <- compute_window_scores(w_df, sessions)
    
    all_window_rows[[length(all_window_rows) + 1]] <- scores %>%
      mutate(
        user_id    = u,
        window_id  = w_name,
        window_idx = as.integer(sub("window_", "", w_name)),
        cohort     = "VAHR",
        start_date = min(w_df$local_date),
        end_date   = max(w_df$local_date)
      )
  }
}

#combine into one df
variability_windows_vahr <- bind_rows(all_window_rows)

#write for external analysis
write.csv(variability_windows_vahr,
          "vahr_window_variability_raw.csv",
          row.names = FALSE)
