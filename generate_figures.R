# Teaching figures for the Results-section workshop (Session 2)
# Kaplan-Meier, ROC and forest plot — clean, chartjunk-free, fully reproducible.
# tidyverse-style ggplot2 + pROC on R 4.4.1; every stochastic step is seeded.

library(survival)
library(ggplot2)
library(pROC)
library(here)

img_dir <- here::here("images")
if (!dir.exists(img_dir)) dir.create(img_dir, recursive = TRUE)

accent  <- "#1f5560"   # deck accent (teal)
accent2 <- "#9e2a2b"   # maroon (.important)
grey    <- "grey55"

theme_talk <- theme_classic(base_size = 15) +
  theme(
    legend.position    = "top",
    legend.title       = element_blank(),
    plot.title.position = "plot",
    plot.subtitle      = element_text(colour = grey, size = 11),
    axis.line          = element_line(colour = "grey40")
  )

# ------------------------------------------------------------------ #
# 1. Kaplan-Meier — two-arm trial, CI bands, censor ticks, log-rank P
# ------------------------------------------------------------------ #
set.seed(2026)
n_arm <- 150
sim_arm <- function(rate, label, n) {
  t_event <- rexp(n, rate = rate)
  t_cens  <- runif(n, 6, 24)                 # administrative censoring (months)
  data.frame(
    time   = pmin(t_event, t_cens),
    status = as.integer(t_event <= t_cens),  # 1 = event, 0 = censored
    arm    = label
  )
}
trial <- rbind(
  sim_arm(0.10, "Treatment", n_arm),
  sim_arm(0.16, "Control",   n_arm)
)
trial$arm <- factor(trial$arm, levels = c("Treatment", "Control"))

fit <- survfit(Surv(time, status) ~ arm, data = trial)

# tidy the survfit into a step data frame (no survminer needed)
arm_names <- sub("arm=", "", names(fit$strata))
km <- data.frame(
  time = fit$time, surv = fit$surv, lower = fit$lower, upper = fit$upper,
  n_censor = fit$n.censor, arm = rep(arm_names, fit$strata)
)
km0 <- data.frame(time = 0, surv = 1, lower = 1, upper = 1, n_censor = 0,
                  arm = arm_names)
km <- rbind(km0, km)
km$arm <- factor(km$arm, levels = c("Treatment", "Control"))
km$lower[is.na(km$lower)] <- km$surv[is.na(km$lower)]
km$upper[is.na(km$upper)] <- km$surv[is.na(km$upper)]

# expand to a proper step ribbon for the CI band
step_ribbon <- function(d) {
  d <- d[order(d$time), ]; n <- nrow(d)
  data.frame(
    time  = as.vector(rbind(d$time[-n],  d$time[-1])),
    lower = as.vector(rbind(d$lower[-n], d$lower[-n])),
    upper = as.vector(rbind(d$upper[-n], d$upper[-n])),
    arm   = d$arm[1]
  )
}
ribbon <- do.call(rbind, lapply(split(km, km$arm), step_ribbon))
ribbon$arm <- factor(ribbon$arm, levels = c("Treatment", "Control"))

sd  <- survdiff(Surv(time, status) ~ arm, data = trial)
p_lr <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p_lab <- if (p_lr < 0.001) "log-rank P < 0.001" else sprintf("log-rank P = %.3f", p_lr)

at <- summary(fit, times = seq(0, 24, 6), extend = TRUE)
risk <- tapply(at$n.risk, sub("arm=", "", at$strata), c)
risk_sub <- sprintf("At risk (mo 0/6/12/18/24)   Treatment: %s    Control: %s",
                    paste(risk[["Treatment"]], collapse = "/"),
                    paste(risk[["Control"]],   collapse = "/"))

p_km <- ggplot() +
  geom_ribbon(data = ribbon, aes(time, ymin = lower, ymax = upper, fill = arm),
              alpha = 0.15) +
  geom_step(data = km, aes(time, surv, colour = arm), linewidth = 1.2) +
  geom_point(data = subset(km, n_censor > 0),
             aes(time, surv, colour = arm), shape = 3, size = 2, stroke = 0.9) +
  annotate("text", x = 1, y = 0.08, hjust = 0, label = p_lab,
           colour = "grey30", size = 4.2) +
  scale_colour_manual(values = c(Treatment = accent, Control = accent2)) +
  scale_fill_manual(values   = c(Treatment = accent, Control = accent2)) +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(c(0, 0.02))) +
  labs(x = "Months since randomisation", y = "Event-free probability",
       title = "Kaplan-Meier estimate by arm", subtitle = risk_sub) +
  theme_talk
ggsave(file.path(img_dir, "km_curve.png"), p_km,
       width = 8.0, height = 5.2, dpi = 200)

# ------------------------------------------------------------------ #
# 2. ROC curve — biomarker vs disease, AUC + Youden-optimal cut-off
# ------------------------------------------------------------------ #
set.seed(7)
n_dx      <- 400
disease   <- rbinom(n_dx, 1, 0.40)
biomarker <- rnorm(n_dx, mean = ifelse(disease == 1, 6.2, 4.8), sd = 1.3)

roc_obj <- roc(disease, biomarker, quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))
best <- coords(roc_obj, "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))

p_roc <- ggroc(roc_obj, legacy.axes = TRUE, colour = accent, linewidth = 1.3) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = grey) +
  annotate("point", x = 1 - best$specificity, y = best$sensitivity,
           colour = accent2, size = 3.4) +
  annotate("text", x = 1 - best$specificity + 0.03, y = best$sensitivity - 0.03,
           hjust = 0, vjust = 1, colour = accent2, size = 4,
           label = sprintf("cut-off = %.1f\nSn %.0f%%, Sp %.0f%%",
                           best$threshold, 100 * best$sensitivity,
                           100 * best$specificity)) +
  annotate("text", x = 0.62, y = 0.12, colour = accent, fontface = "bold",
           size = 5.2, label = sprintf("AUC = %.2f", auc_val)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "1 - Specificity (false-positive rate)",
       y = "Sensitivity (true-positive rate)",
       title = "ROC curve for a continuous biomarker") +
  theme_talk + theme(legend.position = "none")
ggsave(file.path(img_dir, "roc_curve.png"), p_roc,
       width = 6.6, height = 6.2, dpi = 200)

# ------------------------------------------------------------------ #
# 3. Forest plot — subgroup odds ratios with 95% CI (log scale)
# ------------------------------------------------------------------ #
forest_df <- data.frame(
  subgroup = c("Overall", "Age < 50 yr", "Age >= 50 yr",
               "Male", "Female", "Diabetic", "Non-diabetic"),
  or   = c(0.66, 0.71, 0.58, 0.69, 0.62, 0.55, 0.74),
  low  = c(0.48, 0.49, 0.36, 0.47, 0.40, 0.34, 0.50),
  high = c(0.91, 1.03, 0.94, 1.01, 0.96, 0.89, 1.09),
  stringsAsFactors = FALSE
)
forest_df$is_overall <- forest_df$subgroup == "Overall"
forest_df$subgroup <- factor(forest_df$subgroup, levels = rev(forest_df$subgroup))
forest_df$lab <- sprintf("%.2f (%.2f-%.2f)", forest_df$or, forest_df$low, forest_df$high)

p_forest <- ggplot(forest_df, aes(or, subgroup, colour = is_overall)) +
  geom_vline(xintercept = 1, linetype = 2, colour = grey) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.18, linewidth = 0.9) +
  geom_point(aes(shape = is_overall, size = is_overall)) +
  geom_text(aes(x = 1.55, label = lab, fontface = ifelse(is_overall, 2, 1)),
            hjust = 0, size = 4, colour = "grey20") +
  annotate("text", x = 0.55, y = 7.7, label = "favours treatment",
           colour = grey, size = 3.3, hjust = 1) +
  annotate("text", x = 1.05, y = 7.7, label = "favours control",
           colour = grey, size = 3.3, hjust = 0) +
  scale_x_log10(breaks = c(0.3, 0.5, 0.7, 1.0, 1.3),
                limits = c(0.3, 2.6)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.3))) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = c(`FALSE` = accent, `TRUE` = accent2)) +
  scale_shape_manual(values  = c(`FALSE` = 15, `TRUE` = 18)) +
  scale_size_manual(values   = c(`FALSE` = 3, `TRUE` = 4.5)) +
  labs(x = "Adjusted odds ratio (95% CI)", y = NULL,
       title = "Subgroup effects (adjusted OR)") +
  theme_talk +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(colour = "grey92"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(img_dir, "forest_plot.png"), p_forest,
       width = 8.4, height = 5.0, dpi = 200)

cat("Figures written to", img_dir, "\n")
cat(sprintf("  km_curve.png  |  roc_curve.png (AUC %.2f)  |  forest_plot.png\n",
            auc_val))
