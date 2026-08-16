library(dplyr)
library(tidyr)
library(lubridate)
library(circular)
library(purrr)
library(scales)
library(fmsb)
library(ggplot2)

#load files
ios <- read.csv("va_new_ios_logs.csv")
android <- read.csv("va_new_android_logs.csv")

#clean files (with optional filter for VAHR/VAPT)
ios_clean <- ios %>%
  #filter(grepl("^VAPT", association_code)) %>%
  transmute(
    user_id    = users_idx,
    start_time = ymd_hms(start),
    end_time   = ymd_hms(end),
    duration   = duration,
    local_date = as.Date(local_date),
    raw_cat    = category,
    platform   = "ios"
  )

android_clean <- android %>%
  #filter(grepl("^VAPT", association_code)) %>%
  transmute(
    user_id    = users_idx,
    start_time = ymd_hms(start),
    end_time   = ymd_hms(end),
    duration   = duration,
    local_date = as.Date(local_date),
    raw_cat    = genre,
    platform   = "android"
  )

sessions <- bind_rows(ios_clean, android_clean) %>%
  mutate(
    duration   = as.numeric(duration),
    start_hour = hour(start_time) + minute(start_time)/60
  )

sessions <- sessions %>%
  mutate(
    raw_cat = ifelse(raw_cat == "" | is.na(raw_cat), "Unknown", raw_cat)
  )

#specific to the categories seen in these files
category_map <- tribble(
  ~raw_cat,             ~harmonized_category,
  #ios base categories
  "Productivity",      "Productivity",
  "SocialNetworking",  "Social",
  "Finance",           "Finance",
  "PhotoAndVideo",     "PhotoAndVideo",
  "Entertainment",     "Entertainment",
  "Miscellaneous",     "Miscellaneous",
  "Navigation",        "Navigation",
  "Medical",           "Medical",
  "Lifestyle",         "Lifestyle",
  "Weather",           "Weather",
  "Business",          "Business",
  "HealthAndFitness",  "HealthAndFitness",
  "Travel",            "Travel",
  "Utilities",         "Utilities",
  "Games",             "Games",
  "Music",             "Music",
  "FoodAndDrink",      "FoodAndDrink",
  "Sports",            "Sports",
  "News",              "News",
  "Shopping",          "Shopping",
  "Education",         "Education",
  "Books",             "Books",
  "Reference",         "Reference",
  "GraphicsAndDesign", "GraphicsAndDesign",
  "DeveloperTools",    "DeveloperTools",
  
  #android mappings to ios
  "Music & Audio",           "Music",
  "Video Players & Editors", "PhotoAndVideo",
  "Communication",           "Social",
  "Tools",                   "Utilities",
  "Travel & Local",          "Travel",
  "Photography",             "PhotoAndVideo",
  "Health & Fitness",        "HealthAndFitness",
  "Food & Drink",            "FoodAndDrink",
  "Social",                  "Social",
  "News & Magazines",        "News",
  "Books & Reference",       "Books",
  "Maps & Navigation",       "Navigation",
  "Art & Design",            "GraphicsAndDesign",
  "Educational",             "Education",
  
  #android game subtypes
  "Puzzle",        "Games",
  "Role Playing",  "Games",
  "Simulation",    "Games",
  "Card",          "Games",
  "Word",          "Games",
  "Action",        "Games",
  "Casual",        "Games",
  "Strategy",      "Games",
  "Events",        "Games",
  "Board",         "Games",
  "Adventure",     "Games",
  "Casino",        "Games",
  "Trivia",        "Games",
  "Racing",        "Games",
  "Arcade",        "Games",
  
  #any other = miscelaneous
  "Auto & Vehicles",  "Miscellaneous",
  "Personalization",  "Miscellaneous",
  "Beauty",           "Miscellaneous",
  "House & Home",     "Miscellaneous",
  "Dating",           "Miscellaneous",
  "Parenting",        "Miscellaneous",
  "Comics",           "Miscellaneous",
  "Libraries & Demo", "Miscellaneous"
)

sessions <- sessions %>%
  left_join(category_map, by = c("raw_cat" = "raw_cat")) %>%
  mutate(
    harmonized_category = if_else(is.na(harmonized_category),
                                  raw_cat,
                                  harmonized_category)
  )

#shannon entropy helper
shannon_entropy <- function(x) {
  tab <- table(x)
  p   <- tab / sum(tab)
  p   <- p[p > 0]
  -sum(p * log2(p))
}

#split by user and local date
daily <- sessions %>%
  group_by(user_id, local_date) %>%
  summarise(
    total_dur      = sum(duration),
    n_sessions     = n(),
    mean_sess_dur  = mean(duration),
    mean_hour      = mean(start_hour),
    std_hour       = sd(start_hour),
    entropy        = shannon_entropy(harmonized_category),
    .groups = "drop"
  )

#calculate center times for TCV
daily_center <- sessions %>%
  group_by(user_id, local_date) %>%
  summarise(
    center_time = sum(start_hour * duration, na.rm = TRUE) /
      sum(duration, na.rm = TRUE),
    .groups = "drop"
  )

daily <- daily %>%
  left_join(daily_center, by = c("user_id", "local_date"))

#find users with enough days logged
eligible_users <- daily %>%
  count(user_id) %>%
  filter(n >= 21) %>% #can change n to change criteria
  pull(user_id)

daily_eligible <- daily %>%
  filter(user_id %in% eligible_users)

#calculate bvs per user
variability <- daily_eligible %>%
  group_by(user_id) %>%
  summarise(
    VEV_raw = 0.5 * (sd(total_dur)  / mean(total_dur)) +
      0.5 * (sd(n_sessions) / mean(n_sessions)),
    SSV_raw = sd(mean_sess_dur) / mean(mean_sess_dur)
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

variability <- variability %>%
  left_join(TCV_all, by = "user_id")

#ccv helper function: jensen-shannon divergence
jsd <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5*(p + q)

  #avoid dividing by zero
  eps <- 1e-12
  p[p == 0] <- eps
  q[q == 0] <- eps
  m[m == 0] <- eps
  0.5 * sum(p * log2(p/m)) + 0.5 * sum(q * log2(q/m))
}

daily_cat_all <- sessions %>%
  filter(user_id %in% eligible_users) %>%
  group_by(user_id, local_date, harmonized_category) %>%
  summarise(total = sum(duration), .groups = "drop") %>%
  group_by(user_id, local_date) %>%
  mutate(p = total / sum(total)) %>%
  select(user_id, local_date, harmonized_category, p) %>%
  pivot_wider(
    names_from = harmonized_category,
    values_from = p,
    values_fill = 0
  ) %>%
  arrange(user_id, local_date)

CCV_all <- daily_cat_all %>%
  group_by(user_id) %>%
  reframe(
    CCV_raw = {
      if (n() < 2) {
        NA_real_
      } else {
        mats <- as.matrix(cur_data()[, -c(1,2)])
        js_vals <- map_dbl(2:nrow(mats), function(i) {
          jsd(mats[i-1, ], mats[i, ])
        })
        if (is.nan(mean(js_vals))) {
          NA_real_
        } else {
          mean(js_vals)
        }
      }
    },
    .groups = "drop"
  )

### this comment section was a test to change method to cosine divergence
### still working on this...
# jaccard <- function(a, b) {
#   inter <- sum(a & b)
#   union <- sum(a | b)
#   if (union == 0) return(0)
#   1 - inter / union
# }
# 
# CCV_all <- daily_cat_all %>%
#   group_by(user_id) %>%
#   reframe(
#     CCV_raw = {
#       if (n() < 2) {
#         NA_real_
#       } else {
#         mats <- as.matrix(cur_data()[, -c(1,2)])
#         # mats_bin <- (mats > 0) * 1
#         jac_vals <- map_dbl(2:nrow(mats), ~ jaccard(mats[.x-1,], mats[.x]))
#         if (is.nan(mean(jac_vals))) {
#           NA_real_
#         } else {
#           mean(jac_vals)
#         }
#       }
#     },
#     .groups = "drop"
#   )

#add all NaN users to missing_users for reference and single-user plotting
nan_users <- CCV_all %>%
  filter(is.na(CCV_raw)) %>%
  select(user_id)


variability <- variability %>%
  left_join(CCV_all, by = "user_id")

#normalize all to 0-1
variability <- variability %>%
  mutate(
    VEV = percent_rank(VEV_raw),
    TCV = percent_rank(TCV_raw),
    CCV = percent_rank(CCV_raw),
    SSV = percent_rank(SSV_raw),
    BVS = 0.3*VEV + 0.3*TCV + 0.2*CCV + 0.2*SSV
  )

### non-normalized BVS - do not use, just for past reference
# variability <- variability %>%
#   mutate(
#     BVS = 0.3*VEV_raw +
#       0.3*TCV_raw +
#       0.2*CCV_raw +
#       0.2*SSV_raw
#   )


### additional filter for testing whether normalization worked
# odd_users <- variability %>%
#   filter(SSV_raw>1 | TCV_raw>1 | CCV_raw>1 | VEV_raw>1) %>%
#   select(user_id)
# 
# all_users <- variability %>%
#   filter(!is.na(BVS)) %>%
#   pull(user_id)

#calculate BVS for a single user
plottable_users <- variability %>%
  filter(!is.na(BVS)) %>%
  filter(!user_id %in% nan_users$user_id) %>%
  # filter(!user_id %in% odd_users$user_id) %>% #use when checking for normalization
  pull(user_id)

rand <- sample(plottable_users, 10)

#remove for loop and replace i with plottable_users[n] for single user
for (i in rand) {
  scores <- variability %>%
    filter(user_id == i) %>%
    select(VEV, TCV, CCV, SSV)
  
  user_bvs <- variability %>%
    filter(user_id == i) %>%
    pull(BVS)
  
  df <- rbind(
    max = c(1, 1, 1, 1),
    min = c(0, 0, 0, 0),
    scores
  )
  
  message("user ", i, " scores: ", paste0(round(scores, 3), collapse = ", "), " BVS: ", round(user_bvs, 3))
  
  radarchart(df,
             pcol = rgb(0.9, 0.2, 0.5, 0.9), pfcol = rgb(0.9, 0.2, 0.5, 0.5), plwd = 4,
             vlcex = 0.8,
             title = paste("User", i, " - BVS = ", round(user_bvs, 3)),
             sub = paste0("VEV: ", round(scores$VEV, 3),
                           ", TCV: ", round(scores$TCV, 3),
                          ", CCV: ", round(scores$CCV, 3),
                           ", SSV: ", round(scores$SSV, 3))
  )
}

### previous single-user plot code, for reference
# 
# plottable_users <- variability %>%
#   filter(!is.na(BVS)) %>%
#   pull(user_id)
# 
# scores <- variability %>%
#   filter(user_id == plottable_users[5]) %>%
#   select(VEV_raw, TCV_raw, CCV_raw, SSV_raw)
# 
# 
# df <- rbind(
#   max = c(1, 1, 1, 1),
#   min = c(0, 0, 0, 0),
#   scores
# )
# 
# colnames(df) <- c("VEV", "TCV", "CCV", "SSV")
# 
# radarchart(df,
#            pcol = rgb(0.9, 0.2, 0.5, 0.9), pfcol = rgb(0.9, 0.2, 0.5, 0.5), plwd = 4,
#            vlcex = 0.8,
#            title = paste("User", plottable_users[5], "BVS")
# )

### code to save "plottable" scores (non NaN/finite) to a csv for later analysis
# plottable_scores <- variability %>%
#   filter(is.finite(BVS))
# 
# write.csv(plottable_scores, "vahr_scores.csv")

#check for which users are plottable
plottable_users

#functions to calculate and plot every user's BVS
plot_user_bvs <- function(u) {
  
  #get scores
  user_scores <- variability %>%
    filter(user_id == u) %>%
    select(VEV, TCV, CCV, SSV)
  
  user_bvs <- variability %>%
    filter(user_id == u) %>%
    pull(BVS)
  
  #add to df for plotting
  df <- rbind(
    max = c(1, 1, 1, 1),
    min = c(0, 0, 0, 0),
    user_scores
  )
  
  colnames(df) <- c("VEV", "TCV", "CCV", "SSV")
  
  #plot scores
  radarchart(df,
             pcol = rgb(0.9, 0.2, 0.5, 0.9), pfcol = rgb(0.9, 0.2, 0.5, 0.5), plwd = 4,
             vlcex = 0.8,
             title = paste("User", u, " - BVS = ", round(user_bvs, 3)),
             sub = paste0("VEV: ", round(scores$VEV, 3),
                          ", TCV: ", round(scores$TCV, 3),
                          ", CCV: ", round(scores$CCV, 3),
                          ", SSV: ", round(scores$SSV, 3)))
  message(u) #plotting confirmation in R console
}

save_user_bvs_plot <- function(u, out_dir = "working maybe overall abnorm normalized plots") {

  #create output directory if it doesn't exist
  if (!dir.exists(out_dir)) dir.create(out_dir)

  #file path
  file_path <- file.path(out_dir, paste0("user_", u, "_bvs.png"))

  #open PNG device
  png(filename = file_path, width = 800, height = 800)

  #call plotting function
  plot_user_bvs(u)

  #close device
  dev.off()
  #NOTE: PNG device code was entirely Copilot (ln 426-432)
  
  message("Saved plot for user ", u, " → ", file_path) #confirmation in console
}

###loop through and plot every user (uncomment when needed, may take a while)
# for (user in plottable_users) {
#   save_user_bvs_plot(user)
# }

