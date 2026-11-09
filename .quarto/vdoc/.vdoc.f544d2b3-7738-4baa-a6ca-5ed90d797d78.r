#
#
#
#
#
#
#
#
#
#
library(tidyverse)
library(readxl)


# Load the growth data from the Excel file
growth_data <- read.csv("density_growth_data.csv",
    col_types = c("text", "numeric", "text", 
        "text", "text", "numeric", "numeric", 
        "numeric", "numeric", "numeric")) # Specify the column types

survival_data <- read.csv("density_survival_data.csv",
    col_types = c("numeric", "text", "text", 
        "numeric"))

#
#
#

# 1) GROWTH DATA
# A) Convert pixel measurements to mm, and build a UNIQUE container ID.

growth_data <- growth_data %>%
  mutate(
    eye_width_mm     = (eye_width_pixel / pixels) * `known distance` / 1000,
    thorax_width_mm  = (thorax_width_pixel / pixels) * `known distance` / 1000,
    thorax_length_mm = (thorax_length_pixel / pixels) * `known distance` / 1000,
    bin_uid          = interaction(density, bin, drop = TRUE)
  )


# B) Bin-level averages for head width, thorax width, thorax length
bin_averages <- growth_data %>%
  group_by(day, density, bin, bin_uid) %>%
  summarise(
    bin_avg_headwidth     = mean(eye_width_mm, na.rm = TRUE),
    bin_avg_thoraxwidth   = mean(thorax_width_mm, na.rm = TRUE),
    bin_avg_thoraxlength  = mean(thorax_length_mm, na.rm = TRUE),
    .groups = "drop"
  )

# C) Density-level averages (averages of bin averages)
density_averages <- bin_averages %>%
  group_by(day, density) %>%
  summarise(
    density_avg_headwidth    = mean(bin_avg_headwidth, na.rm = TRUE),
    density_avg_thoraxwidth  = mean(bin_avg_thoraxwidth, na.rm = TRUE),
    density_avg_thoraxlength = mean(bin_avg_thoraxlength, na.rm = TRUE),
    .groups = "drop"
  )


# SURVIVAL DATA
#  Remove rows with NA in 'survival'
survival_data <- survival_data %>%
  filter(!is.na(survival))

#  Map density group to initial density, and build the same unique container ID.
initial_density_map <- c(A = 25, B = 50, C = 75, D = 100, E = 125)

survival_data <- survival_data %>%
  mutate(
    initial_density = initial_density_map[density],
    bin_uid         = interaction(density, bin, drop = TRUE)
  )

#  Attach absolute & proportional survival to bin_averages
#    Rename them to bin_abs_survival & bin_prop_survival
bin_averages <- bin_averages %>%
  left_join(
    survival_data %>%
      mutate(
        bin_abs_survival  = survival,
        bin_prop_survival = (survival / initial_density) * 100
      ) %>%
      select(day, density, bin, bin_abs_survival, bin_prop_survival, initial_density),
    by = c("day", "density", "bin")
  )

# Summarize to the density level: average absolute & proportional survival
#    name them density_avg_survival & density_avg_prop_survival
mean_survival_density <- bin_averages %>%
  group_by(day, density) %>%
  summarise(
    density_avg_survival      = mean(bin_abs_survival, na.rm = TRUE),
    density_avg_prop_survival = mean(bin_prop_survival, na.rm = TRUE),
    initial_density           = first(initial_density),
    .groups = "drop"
  )

# Merge these survival summaries into density_averages
density_averages <- density_averages %>%
  left_join(mean_survival_density, by = c("day", "density"))


rm(mean_survival_density)


# CALCULATE YIELD AT THE DENSITY LEVEL

# Add a new column "density_yield" to 'density_averages'
density_averages <- density_averages %>%
  mutate(density_yield = density_avg_survival * density_avg_headwidth^3)   # REVISION: head width CUBED (volumetric proxy)


#
#
#
#
#

# CUSTOM COLORS & LABELS
custom_colours <- c(
  "A" = "#E1D9E1",  # Orange
  "B" = "#CB9E82",  # Green
  "C" = "#AC4C51",  # Blue
  "D" = "#8D2C50",  # Yellow
  "E" = "#3B113B"   # Magenta
)

# Corresponding labels (for legends and facets)
custom_labels <- c(
  "A" = "25",
  "B" = "50",
  "C" = "75",
  "D" = "100",
  "E" = "125"
)
#
#
#
#| label: svg-helper
#| message: false
#| warning: false

# All figures export as SVG into a "figures/" folder via this helper.
library(svglite)

save_svg <- function(plot, filename, width = 6, height = 4) {
  dir.create("figures", showWarnings = FALSE)
  path <- file.path("figures", filename)
  ggsave(path, plot, width = width, height = height,
         device = svglite::svglite, bg = "white")
  cat("Saved:", path, "\n")
}
#
#
#
#
#
#
#
#
#| fig-width: 6
#| fig-height: 4

# Per-day mean + standard error across bins, for each density.
growth_se <- bin_averages %>%
  group_by(density, day) %>%
  summarise(
    mean_hw = mean(bin_avg_headwidth, na.rm = TRUE),
    se_hw   = sd(bin_avg_headwidth, na.rm = TRUE) / sqrt(sum(!is.na(bin_avg_headwidth))),
    .groups = "drop"
  ) %>%
  mutate(lo = mean_hw - se_hw, hi = mean_hw + se_hw)


spline_growth <- growth_se %>%
  group_by(density) %>%
  group_modify(~ {
    sm <- spline(.x$day, .x$mean_hw, n = 200)
    sl <- spline(.x$day, .x$lo,      n = 200)
    sh <- spline(.x$day, .x$hi,      n = 200)
    tibble(day = sm$x, mean_hw = sm$y, lo = sl$y, hi = sh$y)
  }) %>%
  ungroup()

p_growth <- ggplot(density_averages,
                   aes(x = day, y = density_avg_headwidth)) +
  # Error ribbon (±1 SE) around the smoothed mean
  geom_ribbon(data = spline_growth,
              aes(x = day, ymin = lo, ymax = hi,
                  fill = density, group = density),
              alpha = 0.20, colour = NA, inherit.aes = FALSE) +
  # Smooth line through the actual points (spline)
  geom_line(data = spline_growth,
            aes(x = day, y = mean_hw, color = density, group = density),
            linewidth = 0.8, inherit.aes = FALSE) +
  # Points at the real means, black outline
  geom_point(aes(fill = density, group = density),
             shape = 21, size = 2.2, colour = "black", stroke = 0.5) +
  labs(x = "Day", y = "Mean Head Width (mm)") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = custom_colours, labels = custom_labels) +
  scale_fill_manual(values = custom_colours, labels = custom_labels) +
scale_x_continuous(breaks = c(2, 5, 8, 11, 14),
                     limits = c(1, 14),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0.5, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text  = element_text(color = "black")
  )

p_growth
save_svg(p_growth, "growth_by_day_smooth.svg", width = 6, height = 4)
#
#
#
#
#

#| fig-width: 6
#| fig-height: 4


surv_se <- bin_averages %>%
  group_by(density, day) %>%
  summarise(
    mean_s = mean(bin_abs_survival, na.rm = TRUE),
    se_s   = sd(bin_abs_survival, na.rm = TRUE) / sqrt(sum(!is.na(bin_abs_survival))),
    .groups = "drop"
  ) %>%
  mutate(lo = mean_s - se_s, hi = mean_s + se_s)

spline_surv <- surv_se %>%
  group_by(density) %>%
  group_modify(~ {
    sm <- spline(.x$day, .x$mean_s, n = 200)
    sl <- spline(.x$day, .x$lo,     n = 200)
    sh <- spline(.x$day, .x$hi,     n = 200)
    tibble(day = sm$x, mean_s = sm$y, lo = sl$y, hi = sh$y)
  }) %>%
  ungroup()

p_survival <- ggplot(density_averages,
                     aes(x = day, y = density_avg_survival)) +
  geom_ribbon(data = spline_surv,
              aes(x = day, ymin = lo, ymax = hi,
                  fill = density, group = density),
              alpha = 0.20, colour = NA, inherit.aes = FALSE) +
  geom_line(data = spline_surv,
            aes(x = day, y = mean_s, color = density, group = density),
            linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(aes(fill = density, group = density),
             shape = 21, size = 2.2, colour = "black", stroke = 0.5) +
  labs(x = "Day", y = "Mean Survival (Absolute Count)") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = custom_colours, labels = custom_labels) +
  scale_fill_manual(values = custom_colours, labels = custom_labels) +
  scale_x_continuous(breaks = c(2, 5, 8, 11, 14),
                     limits = c(1, 14),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text  = element_text(color = "black")
  )

p_survival


# Proportional survival: same treatment
surv_prop_se <- bin_averages %>%
  group_by(density, day) %>%
  summarise(
    mean_s = mean(bin_prop_survival, na.rm = TRUE),
    se_s   = sd(bin_prop_survival, na.rm = TRUE) / sqrt(sum(!is.na(bin_prop_survival))),
    .groups = "drop"
  ) %>%
  mutate(lo = mean_s - se_s, hi = mean_s + se_s)

spline_surv_prop <- surv_prop_se %>%
  group_by(density) %>%
  group_modify(~ {
    sm <- spline(.x$day, .x$mean_s, n = 200)
    sl <- spline(.x$day, .x$lo,     n = 200)
    sh <- spline(.x$day, .x$hi,     n = 200)
    tibble(day = sm$x, mean_s = sm$y, lo = sl$y, hi = sh$y)
  }) %>%
  ungroup()

p_survival_relative <- ggplot(density_averages,
                              aes(x = day, y = density_avg_prop_survival)) +
  geom_ribbon(data = spline_surv_prop,
              aes(x = day, ymin = lo, ymax = hi,
                  fill = density, group = density),
              alpha = 0.20, colour = NA, inherit.aes = FALSE) +
  geom_line(data = spline_surv_prop,
            aes(x = day, y = mean_s, color = density, group = density),
            linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(aes(fill = density, group = density),
             shape = 21, size = 2.2, colour = "black", stroke = 0.5) +
  labs(x = "Day", y = "Mean Proportional Survival (%)") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = custom_colours, labels = custom_labels) +
  scale_fill_manual(values = custom_colours, labels = custom_labels) +
  scale_x_continuous(breaks = c(2, 5, 8, 11, 14),
                     limits = c(1, 14),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text  = element_text(color = "black")
  )

p_survival_relative

save_svg(p_survival,          "survival_absolute_by_day.svg",     width = 6, height = 4)
save_svg(p_survival_relative, "survival_proportional_by_day.svg", width = 6, height = 4)
#
#
#
#| fig-width: 6
#| fig-height: 4

library(ggpattern)

# Day 1 vs Day 14 yield, using the bin-level yield definition for consistency
# with the line plot (survival × head width^3 per container, averaged).
yield_1_14 <- bin_averages %>%
  mutate(bin_yield = bin_abs_survival * bin_avg_headwidth^3) %>%   # REVISION: head width CUBED (volumetric proxy)
  filter(day %in% c(1, 14)) %>%
  group_by(density, day) %>%
  summarise(mean_yield = mean(bin_yield, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    day       = factor(day, levels = c(1, 14)),
    day_label = if_else(day == "1", "Day 1", "Day 14")
  )

p_yield_1_14 <- ggplot(yield_1_14,
                       aes(x = density, y = mean_yield,
                           fill = density, pattern = day_label)) +
  geom_col_pattern(
    position        = position_dodge(width = 0.8),
    width           = 0.7,
    color           = "black",
    pattern_fill    = "white",
    pattern_color   = "black",
    pattern_density = 0.4
  ) +
  scale_x_discrete(labels = custom_labels) +
  scale_fill_manual(values = custom_colours, guide = "none") +
  scale_pattern_manual(values = c("Day 1" = "stripe", "Day 14" = "none")) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Density Treatment", y = expression("Yield (survivors" %*% "head width"^3*", mm"^3*")")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black"),
    axis.ticks   = element_line(color = "black"),
    axis.text    = element_text(color = "black"),
    legend.title = element_blank()
  )

p_yield_1_14
save_svg(p_yield_1_14, "yield_day1_vs_day14.svg", width = 6, height = 4)


#| fig-width: 6
#| fig-height: 4

# Per-bin yield, then per-day mean + SE across bins, for each density.
yield_se <- bin_averages %>%
  mutate(bin_yield = bin_abs_survival * bin_avg_headwidth^3) %>%   # REVISION: head width CUBED (volumetric proxy)
  group_by(density, day) %>%
  summarise(
    mean_y = mean(bin_yield, na.rm = TRUE),
    se_y   = sd(bin_yield, na.rm = TRUE) / sqrt(sum(!is.na(bin_yield))),
    .groups = "drop"
  ) %>%
  mutate(lo = mean_y - se_y, hi = mean_y + se_y)

spline_yield <- yield_se %>%
  group_by(density) %>%
  group_modify(~ {
    sm <- spline(.x$day, .x$mean_y, n = 200)
    sl <- spline(.x$day, .x$lo,     n = 200)
    sh <- spline(.x$day, .x$hi,     n = 200)
    tibble(day = sm$x, mean_y = sm$y, lo = sl$y, hi = sh$y)
  }) %>%
  ungroup()

p_yield <- ggplot(density_averages,
                  aes(x = day, y = density_yield)) +
  geom_ribbon(data = spline_yield,
              aes(x = day, ymin = lo, ymax = hi,
                  fill = density, group = density),
              alpha = 0.20, colour = NA, inherit.aes = FALSE) +
  geom_line(data = spline_yield,
            aes(x = day, y = mean_y, color = density, group = density),
            linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(aes(fill = density, group = density),
             shape = 21, size = 2.2, colour = "black", stroke = 0.5) +
  labs(x = "Day", y = expression("Yield (survivors" %*% "head width"^3*", mm"^3*")")) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = custom_colours, labels = custom_labels) +
  scale_fill_manual(values = custom_colours, labels = custom_labels) +
  scale_x_continuous(breaks = c(2, 5, 8, 11, 14),
                     limits = c(1, 14),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text  = element_text(color = "black")
  )

p_yield
save_svg(p_yield, "yield_by_day.svg", width = 6, height = 4)
#
#
#
# 4) Stack yield plot  -----------------------------------------------------
library(dplyr)
library(ggplot2)
if (!requireNamespace("ggpattern", quietly = TRUE)) install.packages("ggpattern")
library(ggpattern)

# 1) Filter data to just Day 1 and Day 14 and assign new day labels
yield_1_14 <- density_averages %>%
  filter(day %in% c(1, 14)) %>%
  mutate(
    day = factor(day, levels = c(1, 14)),
    day_label = if_else(day == 1, "Day 1", "Day 14")
  )

# 2) Grouped bar chart with patterns using the updated day labels:
p_yield_1_14 <- ggplot(yield_1_14,
                       aes(x = density,
                           y = density_yield,
                           fill = density,
                           pattern = day_label)) +
  geom_col_pattern(
    position        = position_dodge(width = 0.8),
    width           = 0.7,
    color           = "black",
    pattern_fill    = "white",
    pattern_color   = "black",
    pattern_density = 0.4
  ) +
  scale_x_discrete(labels = custom_labels) +
  scale_fill_manual(values = custom_colours, guide = "none") +
  scale_pattern_manual(
    values = c("Day 1" = "stripe", "Day 14" = "none")
  ) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Density Treatment",
    y = expression("Yield (survivors" %*% "head width"^3*", mm"^3*")")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black"),
    axis.ticks   = element_line(color = "black"),
    axis.text    = element_text(color = "black"),
    legend.title = element_blank()
  )

p_yield_1_14
save_svg(p_yield_1_14, "yield_day1_vs_day14.svg", width = 6, height = 4)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: stats-setup
#| message: false
#| warning: false

library(lme4)       # lmer() / glmer()
library(lmerTest)   # p-values for lmer() via Satterthwaite df
library(emmeans)    # marginal means / trends and Tukey contrasts
library(car)        # Anova()
library(dplyr)
library(tidyr)
library(ggplot2)

ctrl_lmer  <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
ctrl_glmer <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
#
#
#
#
#
#| label: growth-model
#| message: false
#| warning: false

# day_c = day - 7.5 (centred on the middle of the trial). Centring does not
# change the slopes; it just makes the random-slope model fit more stably.
growth_df <- growth_data %>%
  mutate(density = factor(density, levels = c("A", "B", "C", "D", "E")),
         day     = as.numeric(day),
         day_c   = day - 7.5) %>%
  filter(!is.na(eye_width_mm))

# Both fitted with ML so the likelihood-ratio test between them is valid.
growth_ri <- lmer(log(eye_width_mm) ~ density * day_c + (1 | bin_uid),
                  data = growth_df, REML = FALSE, control = ctrl_lmer)
growth_rs <- lmer(log(eye_width_mm) ~ density * day_c + (day_c | bin_uid),
                  data = growth_df, REML = FALSE, control = ctrl_lmer)

cat("\n--- Growth: random-intercept vs random-slope (LRT) ---\n")
print(anova(growth_ri, growth_rs))
cat("Random-slope model singular?", isSingular(growth_rs), "\n")

# Final growth model = random slope (refit with REML for reporting)
growth_model <- update(growth_rs, REML = TRUE)
cat("\n--- Growth: Type III ANOVA (Satterthwaite) ---\n")
print(anova(growth_model))
cat("\n--- Growth: variance components ---\n")
print(VarCorr(growth_model), comp = "Variance")

# Growth slope (log-mm per day) for each density + pairwise contrasts
growth_slopes <- emtrends(growth_model, pairwise ~ density, var = "day_c")
cat("\n--- Growth: slopes per density ---\n")
print(summary(growth_slopes$emtrends, infer = TRUE))
cat("\n--- Growth: pairwise slope differences (Tukey) ---\n")
print(growth_slopes$contrasts)
#
#
#
#
#
#| label: growth-loglinearity
#| fig-width: 7
#| fig-height: 4.5
#| message: false
#| warning: false

# Add a quadratic day term to the SAME model and test whether it improves fit.
growth_quad <- update(growth_rs, . ~ . + I(day_c^2) + density:I(day_c^2))
cat("\n--- Log-linearity: linear vs quadratic (LRT) ---\n")
print(anova(growth_rs, growth_quad))
cat("\n--- Quadratic model ANOVA (does density x day survive?) ---\n")
print(anova(update(growth_quad, REML = TRUE)))

# Visual check: container means of log(head width) with straight line + loess
loglin_bin <- growth_df %>%
  group_by(density, day, bin_uid) %>%
  summarise(log_hw = mean(log(eye_width_mm)), .groups = "drop")
p_loglinearity <- ggplot(loglin_bin, aes(day, log_hw, colour = density)) +
  geom_point(alpha = 0.3, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = custom_colours, labels = custom_labels, name = "Density") +
  labs(x = "Day", y = "log(head width, mm)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"))
p_loglinearity
save_svg(p_loglinearity, "loglinearity_check.svg", width = 7, height = 4.5)

# Residuals vs day: should be flat around 0 if linear on the log scale
growth_df$resid <- residuals(growth_model)
p_resid_day <- ggplot(growth_df, aes(day, resid)) +
  geom_hline(yintercept = 0, colour = "red") +
  geom_point(alpha = 0.1, size = 0.7) +
  geom_smooth(method = "loess", formula = y ~ x, colour = "black") +
  facet_wrap(~ density, nrow = 1, labeller = labeller(density = custom_labels)) +
  labs(x = "Day", y = "Residual (log head width)") +
  theme_minimal(base_size = 12) + theme(panel.grid = element_blank())
p_resid_day
save_svg(p_resid_day, "loglinearity_residuals.svg", width = 8, height = 3)
#
#
#
#
#
#| label: survival-model
#| message: false
#| warning: false

# Cumulative counts are not independent (day 5 contains day 4). So instead we
# model, for each census, how many died SINCE THE PREVIOUS CENSUS out of the
# animals alive at that previous census. Each interval is then a fresh set of
# deaths. Day 0 = the starting number placed in the container.
survival_df <- survival_data %>%
  mutate(density = factor(density, levels = c("A", "B", "C", "D", "E")),
         day = as.numeric(day), alive = survival)

surv_daily <- survival_df %>%
  select(bin_uid, density, bin, initial_density, day, alive) %>%
  bind_rows(distinct(survival_df, bin_uid, density, bin, initial_density) %>%
              mutate(day = 0, alive = initial_density)) %>%
  arrange(bin_uid, day) %>%
  group_by(bin_uid) %>%
  mutate(day_prev = lag(day), alive_prev = lag(alive)) %>%
  ungroup() %>%
  filter(!is.na(alive_prev), alive_prev > 0) %>%
  mutate(interval_days = day - day_prev,
         deaths_raw    = alive_prev - alive,
         deaths        = pmax(deaths_raw, 0),   # a count that went UP is a miscount -> 0 deaths
         survived      = alive_prev - deaths,
         day_c         = day - 7.5)

cat("\n--- Survival intervals ---\n")
cat("Intervals:", nrow(surv_daily), " Containers:", n_distinct(surv_daily$bin_uid), "\n")
cat("Interval lengths (days):\n"); print(table(surv_daily$interval_days))
cat("Intervals where the count went up (set to 0 deaths):", sum(surv_daily$deaths_raw < 0), "\n")

# First try container as a random effect (intercept, then intercept + slope).
# In this daily-mortality data both variance components are estimated at zero
# (singular fits), i.e. once the dependence between successive counts is
# removed, containers do not differ beyond binomial noise. The daily-mortality
# counts are however over-dispersed (ratio ~2.8), so the final model is a
# quasibinomial GLM with no container term. cloglog link + log(interval length)
# offset = discrete-time hazard model: the estimate is the DAILY mortality
# hazard even for the few intervals where a census was missed.
surv_ri <- glmer(cbind(deaths, survived) ~ density * day_c + offset(log(interval_days)) + (1 | bin_uid),
                 data = surv_daily, family = binomial(link = "cloglog"), control = ctrl_glmer)
surv_rs <- glmer(cbind(deaths, survived) ~ density * day_c + offset(log(interval_days)) + (day_c | bin_uid),
                 data = surv_daily, family = binomial(link = "cloglog"), control = ctrl_glmer)
cat("\n--- Survival: container random effects (expect ~0) ---\n")
print(VarCorr(surv_ri)); print(VarCorr(surv_rs))
cat("RI singular?", isSingular(surv_ri), " RS singular?", isSingular(surv_rs), "\n")
cat("Dispersion ratio (RI glmer):",
    round(sum(residuals(surv_ri, type = "pearson")^2) / df.residual(surv_ri), 2), "\n")

# Final model: quasibinomial GLM, F-tests
survival_model <- glm(cbind(deaths, survived) ~ density * day_c + offset(log(interval_days)),
                      data = surv_daily, family = quasibinomial(link = "cloglog"))
cat("\n--- Survival: quasibinomial GLM, dispersion =", round(summary(survival_model)$dispersion, 2), "---\n")
cat("\n--- Survival: Type II F-tests ---\n")
print(car::Anova(survival_model, type = "II", test = "F"))

# Does the density effect change over the trial? (interaction not significant ->
# report density contrasts from the additive model, i.e. averaged over the trial)
survival_add <- glm(cbind(deaths, survived) ~ density + day_c + offset(log(interval_days)),
                    data = surv_daily, family = quasibinomial(link = "cloglog"))
print(car::Anova(survival_add, type = "II", test = "F"))
cat("Day coefficient (log-hazard per day):", round(coef(survival_add)["day_c"], 4),
    " -> hazard multiplies by", round(exp(coef(survival_add)["day_c"]), 3), "per day\n")

# Daily mortality probability per density (mid-trial, per 1 day) + hazard ratios
surv_emm <- emmeans(survival_add, ~ density, type = "response", offset = 0)
cat("\n--- Survival: daily mortality probability per density ---\n")
print(summary(surv_emm))
cat("\n--- Survival: pairwise contrasts (Tukey), cloglog scale. exp(estimate) = hazard ratio ---\n")
surv_pairs <- summary(pairs(surv_emm, adjust = "tukey"))
surv_pairs$hazard_ratio <- exp(-surv_pairs$estimate)   # HR of second-named vs first-named
print(surv_pairs)
#
#
#
#
#
#| label: yield-day14

yield_vol <- bin_averages %>%
  mutate(bin_yield_vol = bin_abs_survival * bin_avg_headwidth^3)

yield_vol_day14 <- yield_vol %>%
  filter(day == 14) %>%
  group_by(density) %>%
  summarise(n_containers = sum(!is.na(bin_yield_vol)),
            mean_yield = mean(bin_yield_vol, na.rm = TRUE),
            se_yield   = sd(bin_yield_vol, na.rm = TRUE) / sqrt(n_containers),
            .groups = "drop") %>%
  mutate(density = custom_labels[as.character(density)]) %>%
  arrange(desc(mean_yield))
cat("\n--- Day-14 volumetric yield proxy (mm^3 units x N), ranked ---\n")
print(yield_vol_day14)

# Same for the old (linear) yield, for comparison
yield_lin_day14 <- bin_averages %>%
  mutate(y = bin_abs_survival * bin_avg_headwidth) %>%
  filter(day == 14) %>% group_by(density) %>%
  summarise(mean_yield = mean(y, na.rm = TRUE), se = sd(y, na.rm=TRUE)/sqrt(sum(!is.na(y))), .groups = "drop") %>%
  arrange(desc(mean_yield))
cat("\n--- (for comparison) Day-14 linear yield proxy ---\n"); print(yield_lin_day14)

# Day 14 proportional survival & head width (descriptive) for the text
cat("\n--- Day-14 descriptives ---\n")
print(bin_averages %>% filter(day == 14) %>% group_by(density) %>%
        summarise(prop_surv = mean(bin_prop_survival, na.rm = TRUE),
                  hw = mean(bin_avg_headwidth, na.rm = TRUE), .groups = "drop"))
#
#
#
growth_df %>% distinct(density, bin) %>% count(density)
#
#
#
#
