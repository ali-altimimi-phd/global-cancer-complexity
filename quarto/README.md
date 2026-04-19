# Quarto Reporting Framework

This directory contains the Quarto source files used by the reporting pipeline.

## Execution Model

Quarto documents in this project are **not intended to be rendered manually**.

All report generation is controlled programmatically via the reporting pipeline:

- `scripts/03_run_report_pipeline.R`
- supporting wrapper functions in `R/wrappers/`

The pipeline determines:
- which templates are used
- which comparisons are rendered
- how outputs are structured

## Directory Contents

This directory includes:

- report templates (`.qmd`)
- supporting Quarto configuration
- legacy and exploratory reporting documents

Not all files are actively used in the current pipeline. Some documents are:

- under development  
- retained for reference  
- part of earlier reporting iterations  

## Rendered Outputs

Rendered reports (HTML, figures, and associated assets) are maintained in a separate repository:

👉 https://github.com/ali-altimimi-phd/cancer-complexity-reports

This separation keeps the codebase lightweight while allowing full browsing of generated reports.

## Notes

- The reporting pipeline should be considered the **authoritative interface** for report generation  
- Direct rendering of individual `.qmd` files may not reflect the intended workflow  