# SQUIRE: Statistical Quality-Assured Integrated Response Estimation

[![CRAN Status](https://www.r-pkg.org/badges/version/SQUIRE)](https://CRAN.R-project.org/package=SQUIRE)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Author:** Richard A. Feiss  
**Version:** 1.0.0  
**License:** MIT  
**Institution:** Minnesota Center for Prion Research and Outreach (MNPRO), University of Minnesota  

---

## Overview

**SQUIRE** (Statistical Quality-Assured Integrated Response Estimation) is an enhanced biological parameter optimization framework that addresses a critical problem in computational biology: **parameter over-interpretation from noisy data**.

Unlike conventional optimizers that attempt parameter fitting on any dataset, SQUIRE implements **statistical gatekeeping** - it first validates whether statistically significant biological effects exist before proceeding with parameter estimation. This prevents computational resources being wasted on noise-fitting and ensures that biological interpretations are statistically justified.

SQUIRE builds upon the established [**GALAHAD**](https://CRAN.R-project.org/package=GALAHAD) optimization framework, adding comprehensive statistical validation and automated parameter type detection.

---

## The Problem SQUIRE Solves

### Traditional Approach (Problematic):
```r
# Standard optimizers fit parameters to ANY data
optimizer(noisy_data) -> parameter_estimates  # Always succeeds!
# Result: False biological interpretations
```

### SQUIRE Approach (Statistically Validated):
```r
# SQUIRE validates BEFORE optimizing
SQUIRE(noisy_data) -> statistical_validation -> {
  if (significant_effects) {
    parameter_optimization -> validated_estimates
  } else {
    "No significant biological effects detected"
  }
}
```

---

## Key Features

###  **Statistical Quality Assurance**
- **Pre-optimization validation**: ANOVA testing for treatment effects
- **Data quality requirements**: Enforced minimum timepoints and replication
- **Effect size assessment**: eta-squared calculation for biological meaningfulness  
- **Negative result handling**: Clear guidance when optimization is unjustified

###  **Biological Intelligence**
- **Automatic parameter typing**: Distinguishes rates, concentrations, bounds
- **Geometry-adaptive optimization**: Applies appropriate methods per parameter type
- **Treatment framework**: Structured comparison against control conditions
- **Response type awareness**: Optimized for germination, growth, survival data

###  **Robust Data Handling**
- **Missing data tolerance**: Adaptive strategies based on data completeness
- **Mixed parameter spaces**: Handles diverse biological model constraints
- **Uncertainty quantification**: Statistical confidence in parameter estimates
- **Reproducible workflows**: Consistent optimization across datasets

---

## Installation

```r
# From CRAN
install.packages("SQUIRE")

# Development version
```

---

## Quick Start

### Basic Usage

```r
library(SQUIRE)

# Load example biological data
data("germination_data")  # Hypothetical dataset

# Statistical quality-assured optimization
results <- SQUIRE(
  data = germination_data,
  treatments = c("Control", "Treatment_A", "Treatment_B"),
  control_treatment = "Control",
  response_type = "germination",
  validation_level = 0.05,
  min_timepoints = 5,
  min_replicates = 3,
  verbose = TRUE
)

# Check results
if (results$optimization_performed) {
  # Significant effects detected - optimization justified
  print(results$parameters)
  print(results$biological_interpretation)
} else {
  # No significant effects - optimization not recommended
  print(results$statistical_advice)
}
```

### Example Output (Positive Result)

```r
# When significant biological effects are detected:
$optimization_performed
[1] TRUE

$statistical_validation
$treatment_effect_pvalue
[1] 0.003

$eta_squared  
[1] 0.74

$parameters
    treatment parameter_1 parameter_2 std_error_1 std_error_2
1     Control      0.12        2.5       0.02        0.3
2 Treatment_A      0.18        3.2       0.03        0.4  
3 Treatment_B      0.24        4.1       0.03        0.5

$biological_interpretation
[1] "Statistically significant treatment effects detected (p=0.003, eta-squared=0.74).
    Parameter optimization justified. Treatment_B shows strongest response."
```

### Example Output (Negative Result)

```r
# When no significant effects are detected:
$optimization_performed
[1] FALSE

$statistical_advice  
[1] "No statistically significant treatment effects detected (p=0.23).
    Consider increasing sample size or re-evaluating experimental design."

$data_quality
$adequate_timepoints: TRUE
$adequate_replication: TRUE
$recommendation: "Insufficient biological signal for parameter optimization"
```

---

## Detailed Workflow

SQUIRE implements a systematic three-stage validation process:

### Stage 1: Data Quality Assessment
-  Minimum timepoints verification (default: 5)
-  Minimum replication verification (default: 3)
-  Data completeness assessment
-  Treatment structure validation

### Stage 2: Statistical Effect Detection  
-  ANOVA for treatment differences
-  Effect size calculation (eta-squared)
-  Statistical significance testing (alpha = 0.05)
-  Biological meaningfulness evaluation

### Stage 3: Validated Parameter Optimization
-  Geometry-adaptive GALAHAD-based optimization
-  Automated parameter type detection
-  Statistical assessment of parameter estimates  
-  Biological interpretation with confidence measures

---

## Supported Response Types

SQUIRE is optimized for biological data patterns:

- **`"germination"`**: Cumulative germination over time
- **`"growth"`**: Plant/organism growth measurements
- **`"survival"`**: Survival analysis with time-to-event data

Each response type uses specialized validation logic and optimization approaches.

---

## Advanced Features

### Custom Validation Levels
```r
# More stringent validation
results <- SQUIRE(
  data = my_data,
  validation_level = 0.01,  # Require p < 0.01
  min_timepoints = 8,       # Require >= 8 timepoints
  min_replicates = 5        # Require >= 5 replicates per treatment
)
```

### Integration with GALAHAD
```r
# Pre-configure GALAHAD parameters (advanced users)
galahad_config <- list(
  geometry_method = "adaptive",
  trust_region_radius = 0.1,
  convergence_tolerance = 1e-6
)

results <- SQUIRE(
  data = my_data,
  galahad_config = galahad_config
)
```

## Citation

When using SQUIRE in publications, please cite:

```
Feiss, R. A. (2025). SQUIRE: Statistical Quality-Assured Integrated Response Estimation. 
R package version 1.0.0. https://CRAN.R-project.org/package=SQUIRE
```

**Please also cite GALAHAD** as SQUIRE builds upon this framework:

```
Feiss, R. A. (2025). GALAHAD: Geometry-Adaptive Lyapunov-Assured Hybrid Optimizer. 
R package version 1.0.0. https://CRAN.R-project.org/package=GALAHAD
```

---

## Development & Support

- **Development**: Minnesota Center for Prion Research and Outreach (MNPRO), University of Minnesota
- **Email**: feiss026@umn.edu

---

## Human-AI Development Transparency

Development followed an **iterative human-machine collaboration**. All algorithmic design, statistical methodologies, and biological validation logic were conceptualized and developed by Richard A. Feiss.

AI systems (*Anthropic Claude*) served as coding and documentation assistants under continuous human oversight, helping with:
- Code optimization and syntax validation
- Statistical method verification  
- Documentation consistency and clarity
- Package compliance checking

AI systems did **not** originate algorithms, statistical approaches, or scientific methodologies.

---

## License

MIT License. See [LICENSE](LICENSE) file for details.