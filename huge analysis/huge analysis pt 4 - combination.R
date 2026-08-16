library(dplyr)
library(ggplot2)

#read overall data csvs produced from previous analysis files
norm_overall <- read.csv("normal_overall_variability_raw.csv")
vahr_overall <- read.csv("vahr_overall_variability_raw.csv")
vapt_overall <- read.csv("vapt_overall_variability_raw.csv")

#combine
overall_all <- bind_rows(norm_overall, vahr_overall, vapt_overall)

#normalize components across all cohorts
overall_all <- overall_all %>%
  mutate(
    VEV = percent_rank(VEV_raw),
    TCV = percent_rank(TCV_raw),
    CCV = percent_rank(CCV_raw),
    SSV = percent_rank(SSV_raw),
    BVS = 0.3 * VEV + 0.3 * TCV + 0.2 * CCV + 0.2 * SSV
  )

overall_all <- overall_all %>%
  filter(is.finite(BVS))
S
write.csv(overall_all, "all_cohorts_overall_normalized.csv", row.names = FALSE)

#read windows
norm_win <- read.csv("normal_window_variability_raw.csv")
vahr_win <- read.csv("vahr_window_variability_raw.csv")
vapt_win <- read.csv("vapt_window_variability_raw.csv")

windows_all <- bind_rows(norm_win, vahr_win, vapt_win)

windows_all <- windows_all %>%
  mutate(
    VEV = percent_rank(VEV_raw),
    TCV = percent_rank(TCV_raw),
    CCV = percent_rank(CCV_raw),
    SSV = percent_rank(SSV_raw),
    BVS = 0.3 * VEV + 0.3 * TCV + 0.2 * CCV + 0.2 * SSV
  )

windows_all <- windows_all %>%
  filter(is.finite(BVS))

write.csv(windows_all, "all_cohorts_windows_normalized.csv", row.names = FALSE)

#boxplots for overall BVS
ggplot(overall_all, aes(x = cohort, y = BVS, fill = cohort)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Overall BVS by Cohort")

ggsave("overall_BVS_boxplot.png", width = 6, height = 4)

#time-series line plots (mean BVS over window_idx) - not properly working at the moment
ggplot(windows_mean, aes(x = window_idx, y = mean_BVS, color = cohort)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  theme_minimal() +
  labs(title = "Mean Sliding-Window BVS Trajectories",
       x = "Window index",
       y = "Mean BVS")

ggsave("windows_BVS_trajectory.png", width = 7, height = 4)

### statistical tests

#ANOVA on overall BVS
anova_res <- aov(BVS ~ cohort, data = overall_all)
summary(anova_res)

#pairwise t-tests (overall BVS)
pairwise_t <- pairwise.t.test(overall_all$BVS, overall_all$cohort)
pairwise_t

#Cohen's d helper function
cohen_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sx <- var(x);    sy <- var(y)
  s_pooled <- sqrt(((nx - 1) * sx + (ny - 1) * sy) / (nx + ny - 2))
  (mean(x) - mean(y)) / s_pooled
}

#effect sizes between cohorts (overall BVS)
bvs_norm <- overall_all$BVS[overall_all$cohort == "normal"]
bvs_vahr <- overall_all$BVS[overall_all$cohort == "VAHR"]
bvs_vapt <- overall_all$BVS[overall_all$cohort == "VAPT"]

d_norm_vahr <- cohen_d(bvs_norm, bvs_vahr)
d_norm_vapt <- cohen_d(bvs_norm, bvs_vapt)
d_vahr_vapt <- cohen_d(bvs_vahr, bvs_vapt)

d_norm_vahr #moderate effect size - noticeable difference
d_norm_vapt #relatively low effect size - not notable difference
d_vahr_vapt #same as norm_vapt