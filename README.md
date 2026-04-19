# 🧬 Cancer Complexity Analysis

### Cancer as a dynamical system: insights from biological complexity and entropy

------------------------------------------------------------------------

## Overview

Cancer is often described as a disease of increasing disorder. This project tests whether that intuition is actually correct.

This repository implements a computational framework for analyzing cancer as a **dynamical system**, where tumorigenesis reflects structured transitions in biological complexity rather than simple randomness.

Rather than treating cancer as uniformly chaotic, the analysis evaluates how **complexity, entropy, and regulatory structure are reorganized** as tissue transitions from normal to malignant states.

Across cancer types, preliminary results suggest a recurring pattern of **regulatory simplification**, accompanied by **functional consolidation and heterogeneous dynamical behavior across biological pathways**, rather than uniform increases in disorder.

------------------------------------------------------------------------

## Tech Stack & Methodology

### Core Environment

-   **Language:** R 4.x (Bioconductor ecosystem)
-   **Workflow:** Configuration-driven, multi-stage analytical pipeline with persistent intermediate artifacts

### Core Analytical Libraries (Bioconductor)

These packages form the backbone of the computational framework:

-   `affy` — RMA normalization and microarray preprocessing\
-   `limma` — differential expression and statistical modeling\
-   `Biobase` — ExpressionSet data structures\
-   `AnnotationDbi` — unified annotation interface

### Data Acquisition & Annotation

-   `GEOquery` — retrieval and parsing of GEO datasets\
-   chip-specific annotation packages (e.g., `hu6800.db`, `hu35ksuba.db`)\
-   GO / KEGG / MSigDB integration via Bioconductor annotation frameworks

### Computational Analysis

-   Custom implementations for:
    -   Shannon entropy
    -   spectral entropy
    -   condition numbers (κ)
    -   effective rank and covariance structure analysis

These are implemented directly in R to allow full control over numerical behavior and interpretation.

### Reporting & Reproducibility

-   **Quarto** — reproducible scientific reporting and document generation\
-   structured logging and configuration-driven execution\
-   modular pipeline design enabling partial recomputation and auditability

### Extensions (In Progress)

-   Python-based latent space modeling:
    -   PyTorch (VAE architectures)
    -   NumPy / pandas for data integration

These components extend the classical framework into learned representations of biological state.

A complete dependency list is available via automated project audit scripts (`R/infrastructure/package_usage.R`).

------------------------------------------------------------------------

## Example Analysis: BLAD/TCC

### Bladder Normal vs Transitional Cell Carcinoma

A representative comparison illustrates the type of structure detected by the framework.

### Summary

| Metric     | Observation                                                  |
|------------|--------------------------------------------------------------|
| Complexity | Net decrease with localized functional gains                 |
| Entropy    | Heterogeneous (coexisting chaotic and anti-chaotic dynamics) |
| Structure  | Selective organization in functional pathways                |

### Interpretation

> In the BLAD/TCC comparison, tumorigenesis reflects a **reorganization of biological complexity**,\
> characterized by loss of regulatory flexibility alongside\
> selective consolidation of functional pathways and heterogeneous dynamical behavior.

This pattern is consistent with a transition toward **structured, task-specific dynamical states**, reflecting regulatory constraint alongside functional specialization.

------------------------------------------------------------------------

### Full Example Output

-   [Example Analysis (BLAD/TCC)](example_blad_tcc.md)
-   [Full BLAD/TCC report](https://ali-altimimi-phd.github.io/cancer-complexity-reports/reports/generated_reports/comparison_report_blad_tcc.html)

This example illustrates the complete analytical workflow, from preprocessing through biological interpretation.

------------------------------------------------------------------------

## Central Question

> **How is biological complexity reorganized during tumorigenesis?**

Rather than assuming monotonic increases in disorder, this project evaluates whether cancer exhibits:

-   **chaotic dynamics** (increased heterogeneity and disorder)
-   **anti-chaotic dynamics** (constraint, clonality, reduced variability)
-   or a **context-dependent interplay between both**

------------------------------------------------------------------------

## Conceptual Framework

### The Cell as a Complex System

Cells exhibit defining properties of complex systems:

-   nonlinear genotype–phenotype relationships
-   network connectivity and feedback
-   sensitivity to initial conditions
-   multiple equilibria and attractor states
-   self-organization and quasi-deterministic behavior

------------------------------------------------------------------------

### Cancer as a Complex Adaptive System

Cancer is modeled as a **state transition process** in gene regulatory space:

-   mutations destabilize regulatory networks
-   systems cross critical thresholds
-   new attractor states emerge (tumor phenotypes)

Importantly:

-   **Solid tumors** often exhibit chaotic deregulation
-   **Hematologic cancers** often exhibit anti-chaotic constraint (clonality)

------------------------------------------------------------------------

### Complexity Loss Hypothesis

A central working hypothesis:

> **Cancer progression is associated with a loss of biological complexity**

This manifests as:

-   reduced effective dimensionality
-   loss of regulatory flexibility
-   convergence toward simplified functional behavior (e.g., proliferation, survival)

------------------------------------------------------------------------

## Data

The analysis is based on the gene expression dataset introduced by Ramaswamy et al. (2001), originally developed for multi-class cancer classification.

-   **Publication:** *Multiclass cancer diagnosis using tumor gene expression signatures* (PNAS, 2001)
-   **GEO accession:** GSE68928
-   **Samples:** \~190 tumor and normal tissue samples
-   **Classes:** 14 distinct tissue and cancer types
-   **Platform:** Affymetrix oligonucleotide microarrays (\~16,000 genes after preprocessing)

Raw **CEL files** are obtained from GEO and processed within this repository, rather than relying on preprocessed expression matrices. This provides full control over normalization, annotation, and filtering.

### Analytical Reframing

Rather than treating the data as a classification problem, the dataset is reorganized into **paired normal vs tumor comparisons** across tissues.

These comparisons are grouped into four major cancer categories:

-   **Carcinomas** (e.g., BLAD/TCC, BR/BRAD, LU/LUAD)
-   **Blastomas** (e.g., Brain/GBM, Brain/MB)
-   **Lymphomas** (e.g., GC/FL, GC/LBCL)
-   **Leukemias** (e.g., PB/AML, PB/B-ALL, PB/T-ALL)

This structure enables:

-   direct comparison of **normal vs malignant state spaces**
-   evaluation of **complexity and entropy changes within tissues**
-   cross-category analysis of **system-level cancer behavior**

### Rationale for Use

The dataset remains well-suited for this analysis because it:

-   spans multiple tissue types within a unified experimental framework
-   provides both **normal and tumor samples** for controlled comparisons
-   captures early evidence of **molecular separability across cancers**

In this project, it serves not as a classification benchmark, but as a **testbed for evaluating hypotheses about biological complexity and dynamical structure in cancer**.

------------------------------------------------------------------------

## Pipeline Architecture

The computational framework is implemented as a **three-stage, configuration-driven analytical system**, designed for reproducibility, modularity, and controlled recomputation.

Each stage operates on well-defined inputs and produces structured outputs that serve as inputs to subsequent stages, forming a **stateful workflow with persistent intermediate artifacts**.

Execution is governed by configuration files, enabling selective activation of stages, parameterized analysis regimes, and consistent reruns across experimental conditions.

------------------------------------------------------------------------

### 1. Preprocessing Pipeline

**Entry point:** `run_preprocessing_pipeline()`

This stage transforms raw microarray data into standardized, analysis-ready representations.

#### Core operations

-   acquisition of raw `.CEL` files (optional, via GEO)
-   chip-specific **RMA normalization**
-   construction of `ExpressionSet` objects (via `Biobase`)
-   integration of GEO metadata (`GEOquery`)
-   metadata cleaning and harmonization:
    -   disease state normalization
    -   tissue label standardization
-   full-chip probe annotation:
    -   Gene Ontology (GO)
    -   KEGG pathways
    -   MSigDB collections

#### Outputs

-   normalized `ExpressionSet` collections (per chip)
-   aligned and cleaned metadata tables
-   probe-to-gene and gene set annotation mappings

This stage establishes the **canonical state space representation** for all downstream analyses.

------------------------------------------------------------------------

### 2. Analysis Pipeline

**Entry point:** `run_analysis_pipeline()`

This stage performs **tissue-specific comparative analysis**, quantifying changes in complexity and entropy between normal and malignant states.

Execution is modular and controlled via configuration flags, supporting multiple analytical regimes.

#### Core components

##### (a) Matrix Construction and Comparison Mapping

-   extraction of tissue-specific expression matrices
-   definition of predefined **normal vs tumor contrasts**
-   alignment of sample groups within each comparison

##### (b) Group-Aware Probe Filtering

-   filtering performed *within comparison context*
-   supported strategies:
    -   **limma-based differential filtering**
    -   **variance-based selection**
-   applied per chip and per comparison

##### (c) Pairwise Comparison Engine

-   execution of all defined tissue comparisons
-   dual analytical modes:
    -   **complexity analysis**
    -   **entropy analysis**
-   evaluation at gene set level:
    -   GO
    -   KEGG
    -   MSigDB

##### (d) Aggregation and Annotation

-   consolidation of gene set–level results
-   attachment of functional annotations
-   harmonization across comparisons and chips

##### (e) Comparison-Level Summarization

-   integration of entropy and complexity outputs
-   computation of derived metrics (e.g., Δ complexity, Δ entropy)
-   generation of unified comparison-level tables

#### Outputs

-   filtered probe sets (per comparison)
-   gene set–level results
-   aggregated datasets
-   comparison-level summary tables

This stage defines the **empirical measurement of dynamical transitions** between biological states.

------------------------------------------------------------------------

### 3. Reporting Pipeline

**Entry point:** `run_report_pipeline()`

This stage converts analytical outputs into structured, interpretable reports.

#### Core operations

-   systematic loading of analysis artifacts
-   generation of complexity and entropy summaries
-   templated construction of comparison-specific narratives
-   GO-based clustering and functional grouping
-   rendering via **Quarto**

#### Outputs

-   structured summary tables (CSV, RDS)
-   comparison-specific reports
-   GO clustering outputs
-   fully rendered Quarto documents

This stage provides the **interpretive layer**, translating quantitative results into biologically meaningful descriptions.

------------------------------------------------------------------------

### Design Characteristics

The pipeline exhibits the following architectural properties:

-   **configuration-driven execution**\
    (explicit control over stages, parameters, and analytical regimes)

-   **stateful workflow design**\
    (intermediate artifacts are persisted and reused)

-   **modular stage separation**\
    (clear boundaries between preprocessing, analysis, and reporting)

-   **support for multiple analytical regimes**\
    (e.g., limma vs variance filtering, entropy vs complexity engines)

-   **fault-tolerant recomputation**\
    (partial pipeline execution without full reruns)

-   **consistent directory and data contracts**\
    (standardized inputs/outputs across stages)

-   **structured logging and traceability**\
    (enabling auditability and reproducibility)

------------------------------------------------------------------------

## Complexity and Entropy Measures

The analysis integrates complementary metrics:

### Complexity (Structure)

-   condition numbers (κ)
-   effective rank
-   sparsity
-   covariance structure

### Entropy (Disorder)

-   Shannon entropy
-   spectral entropy
-   permutation-based significance testing

These enable detection of:

-   complexity loss or gain
-   chaotic vs anti-chaotic transitions
-   system-level reorganization

------------------------------------------------------------------------

## Reports

Rendered reports and supporting HTML outputs are maintained separately from the source code repository.

- **Reports repository:** https://github.com/ali-altimimi-phd/cancer-complexity-reports
- **Published site:** https://ali-altimimi-phd.github.io/cancer-complexity-reports/

This repository contains the source code, pipeline logic, and Quarto reporting framework. The reporting pipeline determines which templates and documents are actively used; rendered outputs are published separately for easier browsing.

------------------------------------------------------------------------

## Project Status

This repository is an **active reconstruction and extension** of prior doctoral research.

-   preprocessing: stable
-   analysis: functional and evolving
-   reporting: expanded beyond original implementation

The repository is presented in curated form for transparency and professional review.

------------------------------------------------------------------------

## Future Directions

This work is being extended using:

-   variational autoencoders (VAE)
-   latent-space geometry
-   modern machine learning approaches

The objective is to connect:

> classical complexity/entropy measures ↔ learned representations of biological state

------------------------------------------------------------------------

## Notes

-   Raw data are not included due to size constraints
-   Scripts assume local data placement or external download
-   Some components are undergoing standardization

------------------------------------------------------------------------

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

------------------------------------------------------------------------

## Citation

If you use this work in academic research, please cite this repository.\
Citation metadata is provided via the `CITATION.cff` file (see “Cite this repository” on GitHub).

------------------------------------------------------------------------

## Acknowledgments

This work extends my doctoral dissertation, *“Chaos and Complexity in Cancer”* (George Mason University, 2004). I am indebted to:

-   **Curtis Jamison**, my dissertation director
-   **Harold J. Morowitz**, whose work on complexity in living systems shaped the conceptual framework of this research
