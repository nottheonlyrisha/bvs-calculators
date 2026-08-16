library(dplyr)
library(tidyr)
library(lubridate)
library(circular)
library(purrr)
library(ggplot2)

### shared category map ###
category_map <- tribble(
  ~raw_cat,             ~harmonized_category,
  # iOS base
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
  
  # Android mappings
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
  
  # Android game subtypes → Games
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
  
  # Other → Miscellaneous
  "Auto & Vehicles",  "Miscellaneous",
  "Personalization",  "Miscellaneous",
  "Beauty",           "Miscellaneous",
  "House & Home",     "Miscellaneous",
  "Dating",           "Miscellaneous",
  "Parenting",        "Miscellaneous",
  "Comics",           "Miscellaneous",
  "Libraries & Demo", "Miscellaneous"
)

### shannon entropy ###
shannon_entropy <- function(x) {
  tab <- table(x)
  p   <- tab / sum(tab)
  p   <- p[p > 0]
  -sum(p * log2(p))
}

### jensen–shannon divergence ###
jsd <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  
  eps <- 1e-12
  p[p == 0] <- eps
  q[q == 0] <- eps
  m[m == 0] <- eps
  
  0.5 * sum(p * log2(p/m)) + 0.5 * sum(q * log2(q/m))
}

### sliding-window helpers ###
get_user_windows <- function(daily_df, user, window_size = 21) {
  user_df <- daily_df %>%
    filter(user_id == user) %>%
    arrange(day_index)
  
  n_days <- nrow(user_df)
  if (n_days < window_size) return(list())
  
  windows <- lapply(1:(n_days - window_size + 1), function(start_idx) {
    user_df %>%
      filter(day_index >= start_idx,
             day_index <  start_idx + window_size)
  })
  names(windows) <- paste0("window_", seq_along(windows))
  windows
}

compute_window_scores <- function(window_df, sessions_df) {
  # VEV + SSV
  vev_ssv <- window_df %>%
    summarise(
      VEV_raw = 0.5 * (sd(total_dur)  / mean(total_dur)) +
        0.5 * (sd(n_sessions) / mean(n_sessions)),
      SSV_raw = sd(mean_sess_dur) / mean(mean_sess_dur)
    )
  
  # TCV
  tcv <- window_df %>%
    summarise(
      TCV_raw = {
        theta <- circular(center_time * 2*pi/24, type = "angles", units = "radians")
        sd.circular(theta)
      }
    )
  
  # CCV
  daily_cat <- sessions_df %>%
    filter(user_id == window_df$user_id[1],
           local_date %in% window_df$local_date) %>%
    group_by(local_date, harmonized_category) %>%
    summarise(total = sum(duration), .groups = "drop") %>%
    group_by(local_date) %>%
    mutate(p = total / sum(total)) %>%
    pivot_wider(names_from = harmonized_category,
                values_from = p,
                values_fill = 0) %>%
    arrange(local_date)
  
  mats <- as.matrix(daily_cat[, -1])
  
  ccv <- tibble(
    CCV_raw = if (nrow(mats) < 2) NA_real_ else {
      js_vals <- map_dbl(2:nrow(mats), function(i) {
        jsd(mats[i-1, ], mats[i, ])
      })
      mean(js_vals)
    }
  )
  
  bind_cols(vev_ssv, tcv, ccv)
}
