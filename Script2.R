# script2.R

# Example: Perform Kaplan-Meier survival analysis
run_survival_analysis <- function(data) {
  library(survival)
  
  # Fit Kaplan-Meier survival curve
  km_fit <- survfit(Surv(OS.time, OS) ~ cancer.type.abbreviation, data = data)
  
  # Print the Kaplan-Meier summary
  print(summary(km_fit))
  
  # Example output
  cat("script2.R: Kaplan-Meier survival analysis completed.\n")
  return(km_fit)
}
