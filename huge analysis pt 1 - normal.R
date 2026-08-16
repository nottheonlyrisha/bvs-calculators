library(dplyr)
library(tidyr)
library(lubridate)
library(circular)
library(purrr)
library(ggplot2)

source("helpers.R")  # or paste shared helpers above

#read files
ios_norm     <- read.csv("va_cohort4_ios_data.csv")
android_norm <- read.csv("va_cohort4_android_data.csv")

#helper function for start and end time formats
parse_min_sec <- function(x) {
  x <- as.character(x)
  ifelse(
    is.na(x) | x == "",
    NA_real_,
    {
      parts <- strsplit(x, ":", fixed = TRUE)[[1]]
      if (length(parts) != 2) {
        NA_real_
      } else {
        mins <- suppressWarnings(as.numeric(parts[1]))
        secs <- suppressWarnings(as.numeric(parts[2]))
        if (is.na(mins) | is.na(secs)) NA_real_ else mins * 60 + secs
      }
    }
  )
}

parse_min_sec_vec <- function(x) {
  vapply(x, parse_min_sec, numeric(1))
}

#clean data
android_clean <- android_norm %>%
  mutate(
    user_id    = users_idx,
    start_sec  = parse_min_sec_vec(start),
    duration   = as.numeric(duration),
    local_date = as.Date(local_date, format = "%d/%m/%Y"),
    raw_cat    = genre,
    platform   = "android",
    start_hour = start_sec / 3600
  ) %>%
  select(user_id, start_hour, duration, local_date, raw_cat, platform)

ios_clean <- ios_norm %>%
  mutate(
    user_id    = users_idx,
    start_sec  = parse_min_sec_vec(start),
    duration   = as.numeric(duration),
    local_date = as.Date(local_date, format = "%d/%m/%Y"),
    raw_cat    = category,
    platform   = "ios",
    start_hour = start_sec / 3600
  ) %>%
  select(user_id, start_hour, duration, local_date, raw_cat, platform)

sessions <- bind_rows(android_clean, ios_clean) %>%
  mutate(raw_cat = ifelse(raw_cat == "" | is.na(raw_cat), "Unknown", raw_cat)) %>%
  left_join(category_map, by = "raw_cat") %>%
  mutate(
    harmonized_category = if_else(is.na(harmonized_category),
                                  raw_cat,
                                  harmonized_category)
  )

#compute daily metrics
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
    center_time = sum(start_hour * duration) / sum(duration),
    .groups = "drop"
  )

daily <- daily %>%
  left_join(daily_center, by = c("user_id", "local_date")) %>%
  group_by(user_id) %>%
  arrange(local_date) %>%
  mutate(day_index = row_number()) %>%
  ungroup()

#compute overall variability for users with over n days of data
eligible_users <- daily %>%
  count(user_id) %>%
  filter(n >= 21) %>% #can change n
  pull(user_id)

daily_eligible <- daily %>%
  filter(user_id %in% eligible_users)

variability_overall_normal <- daily_eligible %>%
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

variability_overall_normal <- variability_overall_normal %>%
  left_join(TCV_all, by = "user_id") %>%
  left_join(CCV_all, by = "user_id") %>%
  mutate(cohort = "normal")

#write data for external analysis
write.csv(variability_overall_normal,
          "normal_overall_variability_raw.csv",
          row.names = FALSE)

#sliding window analysis
window_size <- 21 #can change

all_window_rows <- list()

#loop through every user with enough data and calculate BVS
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
        cohort     = "normal",
        start_date = min(w_df$local_date),
        end_date   = max(w_df$local_date)
      )
  }
}

#combine all calculated BVS windows
variability_windows_normal <- bind_rows(all_window_rows)

#write to external csv for later analysis
write.csv(variability_windows_normal,
          "normal_window_variability_raw.csv",
          row.names = FALSE)
