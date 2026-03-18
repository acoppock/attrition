# Suppress R CMD check notes for ggplot2 NSE column names used in sensitivity_ds
utils::globalVariables(c(
  "p", "upp_est", "low_est", "ci_upper", "ci_lower",
  "value", "label", "hjust", "vjust",
  "change_lower", "change_upper"
))
