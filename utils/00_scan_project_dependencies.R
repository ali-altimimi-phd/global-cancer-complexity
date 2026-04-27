source("R/infrastructure/package_usage.R")

results <- scan_project_packages(project_dir = ".")

print(results$counts, n = Inf)
print(results$counts_by_source, n = Inf)

core_candidates <- identify_core_package_candidates(
  package_counts = results$counts,
  min_count = 5L
)

print(core_candidates, n = Inf)

write_package_usage_outputs(
  scan_results = results,
  output_dir = "output/package_audit",
  prefix = "project"
)

readr::write_csv(
  core_candidates,
  file = "output/package_audit/project_core_package_candidates.csv"
)