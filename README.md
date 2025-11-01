# MICE by mgm — Simulation of Graphical Model–Guided Multiple Imputation for Mixed Data


## 1. Research Overview

This repository extends the methodology introduced in the previous project,  
**“Undergraduate-research-project-1”**, which proposed a  
**Graphical Lasso–guided Multiple Imputation by Chained Equations (Glasso-MICE)** framework.  

While the earlier study focused solely on **continuous data**,  
this project generalizes the framework to **mixed-type data (continuous + categorical)**  
by incorporating the **Mixed Graphical Model (MGM)** into the imputation process.

Specifically, the new framework—**MICE by mgm**—estimates the conditional dependency structure among variables  
using MGM and integrates this structure into the iterative MICE procedure.  
This extension enables a more **realistic and dependency-aware imputation** process for complex datasets  
with heterogeneous variable types.

---

## 2. Research Objectives

- To improve imputation performance under **mixed-data environments**  
- To integrate **MGM-based dependency structures** into the MICE framework  
- To compare the performance of the proposed method with **traditional MICE**  
  and **oracle (complete-data) estimation**, focusing on the accuracy of network recovery  
  measured by **FPR, TPR, and AUC**

---

## 3. Simulation Overview

The simulation is designed to evaluate how accurately each method reconstructs  
the true dependency network (precision matrix) under **10% MCAR missingness** in mixed data.  

| Method | Description |
|--------|--------------|
| **Proposed (MICE by mgm)** | Iterative MICE guided by the network structure estimated via MGM |
| **MICE (Conventional)** | Standard multiple imputation assuming variable independence |
| **Completed Data (Oracle)** | Benchmark estimation using complete data without missingness |

Each method was run across 10 random seeds, and model performance was evaluated using  
**False Positive Rate (FPR)**, **True Positive Rate (TPR)**, and the  
**Area Under the ROC Curve (AUC)** as evaluation metrics.

---

## 4. Workflow

| Step | Description |
|------|--------------|
| ① | Simulation setup and package loading |
| ② | Generation of mixed-type data (continuous + categorical variables) via MGM |
| ③ | Insertion of missing values (MCAR 10%) |
| ④ | Initial imputation using mean/mode substitution |
| ⑤ | Iterative MICE imputation guided by MGM network structure |
| ⑥ | Comparison among Proposed, MICE, and Completed Data methods |
| ⑦ | Calculation of FPR/TPR and ROC curve visualization |

---

## 5. Required Packages

| Package | Purpose |
|----------|----------|
| `mgm` | Estimation of Mixed Graphical Models |
| `mice` | Multiple Imputation by Chained Equations |
| `dplyr` | Data summarization and grouping |
| `caret` | Computation of FPR/TPR via confusion matrices |
| `ggplot2` | Visualization of ROC curves and performance metrics |

---

## 6. Summary of Results

- **Proposed (MICE by mgm)**  
  Incorporating MGM-estimated network structures into the imputation process  
  resulted in **higher TPR and lower FPR**, demonstrating superior structure recovery.  
- **MICE (Conventional)**  
  Assumes variable independence, leading to reduced recovery accuracy.  
- **Completed Data (Oracle)**  
  Serves as the upper bound benchmark for performance comparison.

Performance metrics are summarized as scatter points and mean AUC values  
across multiple simulation runs, with ROC plots generated using `ggplot2`.

---

## 7. How to Run

```r
# Install required packages
install.packages(c("mgm", "mice", "dplyr", "caret", "ggplot2"))

# Run the simulation
source("MICE_by_mgm_Simulation.R")
