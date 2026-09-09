# Suppress R CMD check notes for ggplot2 NSE column names used in sensitivity_ds
utils::globalVariables(c(
  "delta", "estimate_upper", "estimate_lower", "conf.high", "conf.low",
  "value", "label", "hjust", "vjust",
  "change_lower", "change_upper"
))
