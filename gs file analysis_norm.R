library(dplyr)
library(tidyr)
library(lubridate)
library(circular)
library(purrr)
library(scales)
library(fmsb)
library(ggplot2)

#load files
ios <- read.csv("va_cohort4_ios_data.csv")
android <- read.csv("va_cohort4_android_data.csv")

#start and end time format helper
#NOTE: ln 15-37 are entirely Copilot
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

#clean files
android_clean <- android %>%
  mutate(
    user_id    = users_idx,
    start_sec  = parse_min_sec_vec(start),
    end_sec    = parse_min_sec_vec(end),
    duration   = as.numeric(duration),
    local_date = as.Date(local_date, format="%d/%m/%Y"),
    raw_cat    = genre,
    platform   = "android",
    start_hour = start_sec / 3600
  ) %>%
  select(user_id, start_hour, duration, local_date, raw_cat, platform)

ios_clean <- ios %>%
  mutate(
    user_id    = users_idx,
    start_sec  = parse_min_sec_vec(start),
    end_sec    = parse_min_sec_vec(end),
    duration   = as.numeric(duration),
    local_date = as.Date(local_date, format="%d/%m/%Y"),
    raw_cat    = category,
    platform   = "ios",
    start_hour = start_sec / 3600
  ) %>%
  select(user_id, start_hour, duration, local_date, raw_cat, platform)


sessions <- bind_rows(android_clean, ios_clean)


sessions <- sessions %>%
  mutate(raw_cat = ifelse(raw_cat == "" | is.na(raw_cat), "Unknown", raw_cat))

#based on categories in these files
#FUTURE NOTE: figure out how to incorporate user 1020
  #(`genre` is all blank since all apps are in format com.__.android.[app])
  #could maybe exclude individual entirely? however rest of data is viable, only CCV affected
category_map <- tribble(
  ~raw_cat, ~harmonized_category,
  "News", "News",
  "Productivity", "Productivity",
  "Shopping", "Shopping",
  "Utilities", "Utilities",
  "PhotoAndVideo", "PhotoAndVideo",
  "SocialNetworking", "Social",
  "Communication", "Social",
  "Puzzle", "Games",
  "Casino", "Games",
  "Tools", "Utilities",
  "Finance", "Finance",
  "Lifestyle", "Lifestyle",
  "Health & Fitness", "HealthAndFitness",
  "Parenting", "Miscellaneous",
  "Unknown", "Miscellaneous"
)

sessions <- sessions %>%
  left_join(category_map, by=c("raw_cat"="raw_cat")) %>%
  mutate(harmonized_category = ifelse(is.na(harmonized_category),
                                      raw_cat,
                                      harmonized_category))

#shannon entropy helper function
shannon_entropy <- function(x) {
  tab <- table(x)
  p <- tab / sum(tab)
  p <- p[p > 0]
  -sum(p * log2(p))
}

#split by user and local date
daily <- sessions %>%
  group_by(user_id, local_date) %>%
  reframe(
    total_dur     = sum(duration),
    n_sessions    = n(),
    mean_sess_dur = mean(duration),
    mean_hour     = mean(start_hour),
    entropy       = shannon_entropy(harmonized_category)
  )

#compute center times for TCV
daily_center <- sessions %>%
  group_by(user_id, local_date) %>%
  reframe(
    center_time = sum(start_hour * duration) / sum(duration)
  )

daily <- daily %>%
  left_join(daily_center, by=c("user_id","local_date"))

#determine eligible users
eligible_users <- daily %>%
  count(user_id) %>%
  filter(n >= 28) %>% #can change n
  pull(user_id)

daily_eligible <- daily %>%
  filter(user_id %in% eligible_users)

#calculate BVS components
variability <- daily_eligible %>%
  group_by(user_id) %>%
  reframe(
    VEV_raw = 0.5*(sd(total_dur)/mean(total_dur)) +
      0.5*(sd(n_sessions)/mean(n_sessions)),
    SSV_raw = sd(mean_sess_dur)/mean(mean_sess_dur)
  )

TCV_all <- daily_eligible %>%
  group_by(user_id) %>%
  reframe(
    TCV_raw = {
      theta <- circular(center_time * 2*pi/24, type="angles", units="radians")
      sd.circular(theta)
    }
  )

variability <- variability %>%
  left_join(TCV_all, by="user_id")

#ccv jensen-shannon divergence helper function
jsd <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5*(p + q)
  
  eps <- 1e-12
  p[p == 0] <- eps
  q[q == 0] <- eps
  m[m == 0] <- eps
  0.5 * sum(p * log2(p/m)) + 0.5 * sum(q * log2(q/m))
}

daily_cat <- sessions %>%
  filter(user_id %in% eligible_users) %>%
  group_by(user_id, local_date, harmonized_category) %>%
  reframe(total = sum(duration)) %>%
  group_by(user_id, local_date) %>%
  mutate(p = total/sum(total)) %>%
  select(user_id, local_date, harmonized_category, p) %>%
  pivot_wider(names_from=harmonized_category, values_from=p, values_fill=0) %>%
  arrange(user_id, local_date)

CCV_all <- daily_cat %>%
  group_by(user_id) %>%
  reframe(
    CCV_raw = {
      mats <- as.matrix(pick(-c(1,2)))
      if (nrow(mats) < 2) NA_real_ else {
        js_vals <- map_dbl(2:nrow(mats), ~ jsd(mats[.x-1,], mats[.x,]))
        mean(js_vals)
      }
    }
  )

variability <- variability %>%
  left_join(CCV_all, by="user_id")

#normalize all to 0-1
variability <- variability %>%
  mutate(
    VEV = percent_rank(VEV_raw),
    TCV = percent_rank(TCV_raw),
    CCV = percent_rank(CCV_raw),
    SSV = percent_rank(SSV_raw),
    BVS = 0.3*VEV + 0.3*TCV + 0.2*CCV + 0.2*SSV
  )

### non-normalized BVS calculation - ignore, for past reference
# variability <- variability %>%
#   mutate(
#     BVS = 0.3*VEV_raw +
#       0.3*TCV_raw +
#       0.2*CCV_raw +
#       0.2*SSV_raw
#   )

plottable_users <- variability %>%
  filter(!is.na(BVS)) %>%
  pull(user_id)

### write plottable scores to csv for external analysis
# plottable_scores <- variability %>%
#   filter(is.finite(BVS))
# 
# write.csv(plottable_scores, "normal_scores.csv")

#functions to calculate and plot every user's BVS
plot_user_bvs <- function(u) {
  
  #get scores
  user_scores <- variability %>%
    filter(user_id == u) %>%
    select(VEV_raw, TCV_raw, CCV_raw, SSV_raw)
  
  #add to df for plotting
  df <- rbind(
    max = c(1, 1, 1, 1),
    min = c(0, 0, 0, 0),
    user_scores
  )
  
  colnames(df) <- c("VEV", "TCV", "CCV", "SSV")
  
  #plot scores
  radarchart(df,
             pcol = rgb(0.2, 0.5, 0.5, 0.9), pfcol = rgb(0.2, 0.5, 0.5, 0.5), plwd = 4,
             vlcex = 0.8,
             title = paste("User", u, "BVS")
  )
  print(u)
  print(user_scores) #confirmation in console
}

save_user_bvs_plot <- function(u, out_dir = "normal polar plots (28 days)") {

  #create output directory if it doesn't exist
  if (!dir.exists(out_dir)) dir.create(out_dir)

  #file path
  file_path <- file.path(out_dir, paste0("user_", u, "_bvs.png"))

  #open PNG device
  png(filename = file_path, width = 800, height = 800)

  #call your existing plotting function
  plot_user_bvs(u)

  #close device
  dev.off()
  #NOTE: PNG device code was entirely Copilot (ln 261-270)

  message("Saved plot for user ", u, " → ", file_path) #confirmation in console
}

###loop through and plot every user (uncomment when needed)
# for (user in plottable_users) {
#   plot_user_bvs(user)
# }

