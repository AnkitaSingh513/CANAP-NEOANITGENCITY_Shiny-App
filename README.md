## CANAP: Cancer Neo-Antigen Prediction Shiny App

## Overview

**CANAP (Cancer Neo-Antigen Prediction)** is a Shiny-based application developed to aid cancer research by analyzing and visualizing mutation data. It empowers researchers with interactive tools to filter and explore gene mutations, visualize mutation frequencies, perform Gene Ontology (GO) analysis, and compare mutation patterns across multiple cancer types.

This platform supports research in cancer neo-antigenicity, providing actionable insights for precision oncology and immunotherapy.

---

## Key Features

- **Interactive Filtering**: Dynamically filter datasets by genes, cancer types, mutation thresholds, and GO terms.
- **Visualization Tools**:
  - Bar plots and pie charts for mutation distributions.
  - Scatter plots for mutation frequency vs. log difference analysis.
  - Heatmaps for mutation patterns across top genes.
- **GO Analysis**:
  - Explore enriched GO terms for genes with significant mutations.
  - Visualize gene-GO relationships with interactive heatmaps.
- **Comparative Analysis**: Compare top mutated genes across cancer types to identify shared or unique patterns.
- **Data Export**: Download filtered datasets for downstream analysis.

---

## Getting Started

### Prerequisites

Ensure you have the following installed on your system:
- R (version 4.0 or later)
- RStudio (optional but recommended)
- Required R packages:
  - `shiny`
  - `shinydashboard`
  - `plotly`
  - `data.table`
  - `dplyr`
  - `tidyr`
  - `reshape2`
  - `DT`
  - `RColorBrewer`
  - `shinycssloaders`

To install all required packages, run:
```R
install.packages(c("shiny", "shinydashboard", "plotly", "data.table", "dplyr", 
                   "tidyr", "reshape2", "DT", "RColorBrewer", "shinycssloaders"))
