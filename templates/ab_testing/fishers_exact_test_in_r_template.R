
# Fisher's Exact Test in R

# Summary:
# Comparing two landing pages: control vs. treatment.
# Original N = 30,000, but truncated to 305 users.
# Fisher's test used due to small counts.
# Treatment shows statistically significant improvement (α = 0.05), small practical effect (Cohen’s h = 0.02).

# ----------------------------
# 📦 Setup
# ----------------------------

# Install packages
install.packages("exact2x2")
install.packages("statmod")

# Load libraries
library(exact2x2)
library(statmod)
library(glue)

# ----------------------------
# 🧪 Test Design – Parameters & Hypothesis
# ----------------------------

alpha <- 0.05
power <- 0.80
control <- 0.14
effect <- 0.05
mde <- control * effect
treatment <- control + mde

print(paste('Control:', control))
print(paste('Treatment:', treatment))

p_1 <- treatment
p_2 <- control
p1_label <- "Treatment"
p2_label <- "Control"

alternative <- "greater"

hypothesis <- switch(alternative,
  greater = sprintf("%s (%.4f) is greater than %s (%.4f)", p1_label, p_1, p2_label, p_2),
  less = sprintf("%s (%.4f) is less than %s (%.4f)", p1_label, p_1, p2_label, p_2),
  two.sided = sprintf("%s (%.4f) is different from %s (%.4f)", p1_label, p_1, p2_label, p_2)
)

cat("Hypothesis:", hypothesis)

# ----------------------------
# 📏 Effect Size – Cohen’s h
# ----------------------------

proportion_effectsize <- function(treatment, control) {
  2 * asin(sqrt(treatment)) - 2 * asin(sqrt(control))
}

effect_size <- proportion_effectsize(treatment, control)

cat(sprintf("Control = %.4f
", control))
cat(sprintf("Treatment = %.4f
", treatment))
cat(sprintf("Minimum Detectable Effect (MDE): %.3f
", mde))
cat(sprintf("Effect Size (Cohen's h): %.3f
", effect_size))

# ----------------------------
# 🧮 Sample Size Simulation and Finder
# ----------------------------

simulate_fisher_power <- function(p1, p2, n1, n2, alpha, reps = 1000, alternative, seed=100) {
  set.seed(seed)
  rejects <- replicate(reps, {
    x1 <- rbinom(1, n1, p1)
    x2 <- rbinom(1, n2, p2)
    tbl <- matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2, byrow = TRUE)
    fisher.test(tbl, alternative = alternative)$p.value < alpha
  })
  mean(rejects)
}

n_1 <- 30000
n_2 <- 30000

estimated_power <- simulate_fisher_power(p1=p_1, p2=p_2, n_1, n_2, alpha=alpha, alternative=alternative)
cat(sprintf("Power for manual estimate: %.3f
", estimated_power))

find_min_sample_size <- function(p1, p2, alpha, power, max_n = 50000,
                                 reps = 1000, alternative, step = 1000, seed=100) {
  set.seed(seed)
  for (n in seq(1000, max_n, by = step)) {
    sim_power <- simulate_fisher_power(p1, p2, n, n, alpha, reps, alternative, seed)
    cat(sprintf("n = %d → power = %.3f
", n, sim_power))
    if (sim_power >= power) return(n)
  }
  return(NA)
}

set.seed(100)
cat("
Searching for minimum required sample size...
")
min_n <- find_min_sample_size(p1 = p_1, p2 = p_2, alpha = alpha, power = power, alternative = alternative)
cat(sprintf("
Minimum sample size per group to achieve %.0f%% power: %d
", power * 100, min_n))

# ----------------------------
# 📥 Input Data
# ----------------------------

control_conversions <- 7
treatment_conversions <- 18
control_no_conversions <- 150
treatment_no_conversions <- 130

print(p1_label)
print(p2_label)

# ----------------------------
# 📊 Contingency Table
# ----------------------------

table <- matrix(c(control_conversions, control_no_conversions, treatment_conversions, treatment_no_conversions), nrow = 2, byrow = TRUE)
colnames(table) <- c("Converted", "Not_Converted")
rownames(table) <- c("Control", "Treatment")
print(table)

if (rownames(table)[1] != p1_label) {
  table_to_use <- table[c(2, 1), ]
} else {
  table_to_use <- table
}
print(table_to_use)

# ----------------------------
# 📈 Conversion Rates
# ----------------------------

n1 <- sum(table_to_use[1, ])
n2 <- sum(table_to_use[2, ])

p1 <- table_to_use[1, "Converted"] / n1
p2 <- table_to_use[2, "Converted"] / n2

groups <- c(p1 = p1_label, p2 = p2_label)

print(glue("p1: ", "{groups['p1']} Conversion Rate: {round(p1 * 100, 2)}%"))
print(glue("p2: ", "{groups['p2']} Conversion Rate: {round(p2 * 100, 2)}%"))

result_hypothesis <- switch(alternative,
  greater = sprintf("%s (%.4f) is greater than %s (%.4f)", p1_label, p1, p2_label, p2),
  less = sprintf("%s (%.4f) is less than %s (%.4f)", p1_label, p1, p2_label, p2),
  two.sided = sprintf("%s (%.4f) is different from %s (%.4f)", p1_label, p1, p2_label, p2)
)

cat("Result Hypothesis:", result_hypothesis)

# ----------------------------
# 📏 Effect Size (from observed data)
# ----------------------------

abs_diff <- abs(p1 - p2)

proportion_effectsize <- function(control, treatment) {
  2 * asin(sqrt(treatment)) - 2 * asin(sqrt(control))
}

h <- proportion_effectsize(control, treatment)

cat(sprintf("Absolute difference: %.3f (%.1f%%)\n", abs_diff, abs_diff * 100))
cat(sprintf("Cohen's h: %.3f\n", h))

interpret_h <- function(h) {
  if (abs(h) < 0.2) return("negligible")
  if (abs(h) < 0.5) return("small")
  if (abs(h) < 0.8) return("medium")
  return("large")
}
cat(sprintf("Effect size interpretation: %s\n", interpret_h(h)))

# ----------------------------
# 🧪 Fisher’s Exact Test
# ----------------------------

print(table_to_use)
print(alternative)
print(alpha)
print(power)

result <- exact2x2(table_to_use, alternative = alternative, conf.level = 1 - alpha, tsmethod = "central")
print(result)

p_value <- result$p.value
print(paste("p-value:", round(p_value, 3)))

pvalue_message <- if (p_value < alpha) {
  sprintf("Because the p-value (%.3f) is less than alpha (%.3f), this result is statistically significant at the %.0f%% confidence level.",
          p_value, alpha, (1 - alpha) * 100)
} else {
  sprintf("Because the p-value (%.3f) is greater than or equal to alpha (%.3f), this result is not statistically significant at the %.0f%% confidence level.",
          p_value, alpha, (1 - alpha) * 100)
}

cat(strwrap(pvalue_message, width = 80), sep = "\n")

# ----------------------------
# 🎯 Confidence Interval Interpretation
# ----------------------------

lower_ci <- result$conf.int[1]
upper_ci <- result$conf.int[2]

if (alternative == "less") {
  if (upper_ci < 1) {
    percent_diff <- (1 - upper_ci) * 100
    sentence <- paste0(
      sprintf("95%% CI Upper Bound for Odds Ratio: %.3f.", upper_ci), "\n\n",
      sprintf("With %.0f%% confidence, the control group (p₁) has lower odds of conversion than the treatment group (p₂).", (1 - alpha) * 100), "\n",
      sprintf("The odds of conversion in the control group are up to %.1f%% lower than in the treatment group.", percent_diff), "\n",
      "This supports the hypothesis that treatment is better than control."
    )
  } else {
    sentence <- paste0(
      sprintf("95%% CI Upper Bound for Odds Ratio: %.3f.", upper_ci), "\n\n",
      sprintf("With %.0f%% confidence, we cannot rule out that the treatment group is not better than the control group.", (1 - alpha) * 100), "\n",
      "This does not support the hypothesis that treatment is better."
    )
  }

} else if (alternative == "greater") {
  if (lower_ci > 1) {
    percent_diff <- (lower_ci - 1) * 100
    sentence <- paste0(
      sprintf("95%% CI Lower Bound for Odds Ratio: %.3f.", lower_ci), "\n\n",
      sprintf("With %.0f%% confidence, the treatment group (p₂) has higher odds of conversion than the control group (p₁).", (1 - alpha) * 100), "\n",
      sprintf("The odds of conversion in the treatment group are at least %.1f%% higher than in the control group.", percent_diff), "\n",
      "This supports the hypothesis that treatment is better than control."
    )
  } else {
    sentence <- paste0(
      sprintf("95%% CI Lower Bound for Odds Ratio: %.3f.", lower_ci), "\n\n",
      sprintf("With %.0f%% confidence, we cannot rule out that the treatment group is not better than the control group.", (1 - alpha) * 100), "\n",
      "This does not support the hypothesis that treatment is better."
    )
  }

} else if (alternative == "two.sided") {
  if (lower_ci > 1) {
    percent_diff <- (lower_ci - 1) * 100
    sentence <- paste0(
      sprintf("95%% CI: [%.3f, %.3f].", lower_ci, upper_ci), "\n\n",
      sprintf("With %.0f%% confidence, the treatment group (p₂) has higher odds of conversion than the control group (p₁).", (1 - alpha) * 100), "\n",
      sprintf("The odds of conversion in the treatment group are at least %.1f%% higher than in the control group.", percent_diff), "\n",
      "This supports a significant difference favoring treatment."
    )
  } else if (upper_ci < 1) {
    percent_diff <- (1 - upper_ci) * 100
    sentence <- paste0(
      sprintf("95%% CI: [%.3f, %.3f].", lower_ci, upper_ci), "\n\n",
      sprintf("With %.0f%% confidence, the control group (p₁) has lower odds of conversion than the treatment group (p₂).", (1 - alpha) * 100), "\n",
      sprintf("The odds of conversion in the control group are up to %.1f%% lower than in the treatment group.", percent_diff), "\n",
      "This supports a significant difference favoring treatment."
    )
  } else {
    sentence <- paste0(
      sprintf("95%% CI: [%.3f, %.3f].", lower_ci, upper_ci), "\n\n",
      sprintf("With %.0f%% confidence, we cannot rule out no difference in odds between treatment and control.", (1 - alpha) * 100), "\n",
      "This does not support a statistically significant difference."
    )
  }
}

cat(sentence, "\n\n")

includes_one <- (lower_ci <= 1) && (upper_ci >= 1)

message <- if (includes_one) {
  sprintf("Because the interval includes 1, this result is not statistically significant at the %.0f%% confidence level.", (1 - alpha) * 100)
} else {
  sprintf("Because the interval does not include 1, this result is statistically significant at the %.0f%% confidence level.", (1 - alpha) * 100)
}

cat(message)

# ----------------------------
# ⚡ Statistical Power of the Result
# ----------------------------

print(table_to_use)

set.seed(100)
result_power <- power.fisher.test(n1 = n1, n2 = n2, p1 = p1, p2 = p2,
                                  alpha = alpha,
                                  alternative = alternative,
                                  nsim = 10000)

power_pct <- sprintf("%.1f", result_power * 100)
cat("Result Power:", power_pct, "%\n\n")

power_sentence <- if (result_power < 0.8) {
  sprintf("Our test was underpowered (~%s%% power). There was a higher chance we failed to detect a true difference.", power_pct)
} else {
  sprintf("Our test was adequately powered (~%s%% power).", power_pct)
}
cat(power_sentence, "\n")

# ----------------------------
# 🧾 Results Summary
# ----------------------------

print(table_to_use)

print(glue("p1: ","{groups['p1']} Conversion Rate: {round(p1 * 100, 2)}%"))
print(glue("p2: ","{groups['p2']} Conversion Rate: {round(p2 * 100, 2)}%"))

cat("Result Hypothesis:", result_hypothesis, "\n")

cat(sprintf("Absolute difference: %.3f (%.1f%%)\n", abs_diff, abs_diff * 100))
cat(sprintf("Cohen's h: %.3f\n", h))
cat(sprintf("Effect size interpretation: %s\n", interpret_h(h)))

print(result)
cat(strwrap(pvalue_message, width = 80), sep = "\n")

cat(sentence, "\n")
cat(message, "\n")

cat("Result Power:", power_pct, "%\n\n")
cat(power_sentence, "\n")

# ----------------------------
# 🧾 Odds Ratio Conversion
# ----------------------------

print(paste("Control:",control))
print(paste("Treatment:",treatment))

control_odds=control/(1-control)
treatment_odds=treatment/(1-treatment)
odds_ratio=treatment_odds/control_odds

print(paste("Odds Ratio:",round(odds_ratio,2)))
print(paste("Lower CI:",round(lower_ci,2)))
# print(paste("Upper CI:",round(upper_ci,2)))

# Your hypothesized Odds Ratio (OR) of 1.06 (based on Relative Rate (RR) = 0.147/0.14 =1.05)
# 
# Control odds = 0.14 / (1 - 0.14) = 0.1628
# Treatment odds = 0.147 / (1 - 0.147) = 0.1724
# Treatment odds / Control odds = 0.1724 / 0.1628 = 1.059
# 5% increase in probability leads to a 6% increase in odds
# 
# If your minimum meaningful effect was an odds ratio of 1.06, then:
#   
#   You’d want your lower CI bound ≥ 1.06
# 
# Your result of 1.295 >1.06 → it's statistically significant, and stronger than hypothesized
