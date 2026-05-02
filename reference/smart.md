# Imputed SMART study data

Imputed SMART study data (no missing values).

## Usage

``` r
data(smart)
```

## Format

A numeric matrix with 3873 samples, and 29 variables (27 variables +
time variable + event variable):

- Demographics

  - SEX - gender

  - AGE - age in years

- Classical risk factors

  - SMOKING - smoking (never, former, current)

  - PACKYRS - in years

  - ALCOHOL - alcohol use (never, former, current)

  - BMI - Body mass index, in kg/m^2

  - DIABETES

- Blood pressure

  - SYSTH - Systolic, by hand, in mm Hg

  - SYSTBP - Systolic, automatic, in mm Hg

  - DIASTH - Diastolic, by hand, in mm Hg

  - DIASTBP - Diastolic, automatic, in mm Hg

- Lipid levels

  - CHOL - Total cholesterol, in mmol/L

  - HDL - High-density lipoprotein cholesterol, in mmol/L

  - LDL - Low-density lipoprotein cholesterol, in mmol/L

  - TRIG - Triglycerides, in mmol/L

- Previous symptomatic atherosclerosis

  - CEREBRAL - Cerebral

  - CARDIAC - Coronary

  - PERIPH - Peripheral

  - AAA - Abdominal aortic aneurysm

- Markers of atherosclerosis

  - HOMOC - Homocysteine, in \\\mu\\mol/L

  - GLUT - Glutamine, in \\\mu\\mol/L

  - CREAT - Creatinine clearance, in mL/min

  - ALBUMIN - Albumin (no, micro, macro)

  - IMT - Intima media thickness, in mm

  - STENOSIS - Carotid artery stenosis \> 50%

## Note

See `data-raw/smart.R` for the code to generate this data.

## References

Steyerberg, E. W. (2008). Clinical prediction models: a practical
approach to development, validation, and updating. Springer Science &
Business Media.

## Examples

``` r
data(smart)
dim(smart)
#> [1] 3873   29
```
