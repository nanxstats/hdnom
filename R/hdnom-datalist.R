#' Original SMART Study Data
#'
#' Original SMART study data (with missing values) from Steyerberg et, al. 2008.
#'
#' @docType data
#' @name smarto
#' @usage data(smarto)
#'
#' @format
#' A numeric matrix with 3873 samples, and 29 variables
#' (27 variables + time variable + event variable):
#' \itemize{
#' \item Demographics
#' \itemize{
#' \item SEX - gender
#' \item AGE - age in years
#' }
#' \item Classical risk factors
#' \itemize{
#' \item SMOKING - smoking (never, former, current)
#' \item PACKYRS - in years
#' \item ALCOHOL - alcohol use (never, former, current)
#' \item BMI - Body mass index, in kg/m^2
#' \item DIABETES
#' }
#' \item Blood pressure
#' \itemize{
#' \item SYSTH - Systolic, by hand, in mm Hg
#' \item SYSTBP - Systolic, automatic, in mm Hg
#' \item DIASTH - Diastolic, by hand, in mm Hg
#' \item DIASTBP - Diastolic, automatic, in mm Hg
#' }
#' \item Lipid levels
#' \itemize{
#' \item CHOL - Total cholesterol, in mmol/L
#' \item HDL - High-density lipoprotein cholesterol, in mmol/L
#' \item LDL - Low-density lipoprotein cholesterol, in mmol/L
#' \item TRIG - Triglycerides, in mmol/L
#' }
#' \item Previous symptomatic atherosclerosis
#' \itemize{
#' \item CEREBRAL - Cerebral
#' \item CARDIAC - Coronary
#' \item PERIPH - Peripheral
#' \item AAA - Abdominal aortic aneurysm
#' }
#' \item Markers of atherosclerosis
#' \itemize{
#' \item HOMOC - Homocysteine, in \eqn{\mu}{u}mol/L
#' \item GLUT - Glutamine, in \eqn{\mu}{u}mol/L
#' \item CREAT - Creatinine clearance, in mL/min
#' \item ALBUMIN - Albumin (no, micro, macro)
#' \item IMT - Intima media thickness, in mm
#' \item STENOSIS - Carotid artery stenosis > 50\%
#' }
#' }
#'
#' @references
#' Steyerberg, E. W. (2008). Clinical prediction models:
#' a practical approach to development, validation, and updating.
#' Springer Science & Business Media.
#'
#' @examples
#' data(smarto)
#' dim(smarto)
NULL

#' Imputed SMART Study Data
#'
#' Imputed SMART study data (no missing values).
#'
#' @docType data
#' @name smart
#' @usage data(smart)
#'
#' @format
#' A numeric matrix with 3873 samples, and 29 variables
#' (27 variables + time variable + event variable):
#' \itemize{
#' \item Demographics
#' \itemize{
#' \item SEX - gender
#' \item AGE - age in years
#' }
#' \item Classical risk factors
#' \itemize{
#' \item SMOKING - smoking (never, former, current)
#' \item PACKYRS - in years
#' \item ALCOHOL - alcohol use (never, former, current)
#' \item BMI - Body mass index, in kg/m^2
#' \item DIABETES
#' }
#' \item Blood pressure
#' \itemize{
#' \item SYSTH - Systolic, by hand, in mm Hg
#' \item SYSTBP - Systolic, automatic, in mm Hg
#' \item DIASTH - Diastolic, by hand, in mm Hg
#' \item DIASTBP - Diastolic, automatic, in mm Hg
#' }
#' \item Lipid levels
#' \itemize{
#' \item CHOL - Total cholesterol, in mmol/L
#' \item HDL - High-density lipoprotein cholesterol, in mmol/L
#' \item LDL - Low-density lipoprotein cholesterol, in mmol/L
#' \item TRIG - Triglycerides, in mmol/L
#' }
#' \item Previous symptomatic atherosclerosis
#' \itemize{
#' \item CEREBRAL - Cerebral
#' \item CARDIAC - Coronary
#' \item PERIPH - Peripheral
#' \item AAA - Abdominal aortic aneurysm
#' }
#' \item Markers of atherosclerosis
#' \itemize{
#' \item HOMOC - Homocysteine, in \eqn{\mu}{u}mol/L
#' \item GLUT - Glutamine, in \eqn{\mu}{u}mol/L
#' \item CREAT - Creatinine clearance, in mL/min
#' \item ALBUMIN - Albumin (no, micro, macro)
#' \item IMT - Intima media thickness, in mm
#' \item STENOSIS - Carotid artery stenosis > 50\%
#' }
#' }
#'
#' @references
#' Steyerberg, E. W. (2008). Clinical prediction models:
#' a practical approach to development, validation, and updating.
#' Springer Science & Business Media.
#'
#' @examples
#' data(smart)
#' dim(smart)
#'
#' # # Code for generating the smart dataset by imputing the smarto data
#' # library("mice")
#' # library("Hmisc")
#' # data(smarto)
#' # set.seed(1000)
#' # smarti = aregImpute(~ I(TEVENT) + EVENT + SEX + I(AGE) +
#' #   DIABETES + CEREBRAL + CARDIAC + AAA + PERIPH +
#' #   STENOSIS + SYSTBP + DIASTBP + SYSTH + DIASTH +
#' #   I(LENGTH) + I(WEIGHT) + I(BMI) +
#' #   I(CHOL) + I(HDL) + I(LDL) + I(TRIG) + I(HOMOC) + I(GLUT) + I(CREAT) + I(IMT) +
#' #   as.factor(ALBUMIN) + as.factor(SMOKING) + I(PACKYRS) + as.factor(ALCOHOL),
#' #   n.impute = 10, data = smarto)
#' # imputed = impute.transcan(smarti, imputation = 1, data = smarto,
#' #   list.out = TRUE, pr = FALSE, check = FALSE)
#' # smart = smarto
#' # smart[names(imputed)] = imputed
#' # str(smart)
NULL
