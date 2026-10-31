#' SIGA: Subgroup Identification Based on a Genetic Algorithm
#'
#' @description
#' The SIGA method searches for a subgroup of patients, defined by a
#' conjunction of biomarker conditions, in which the treatment effect is
#' considerably larger than in the complementary population. It is intended for
#' exploratory analyses of two-arm randomised clinical trials with
#' time-to-event, binary or continuous endpoints.
#'
#' The procedure (Zhao et al.) consists of four steps:
#' \enumerate{
#'   \item **Prescreening** of candidate biomarkers by the treatment-by-biomarker
#'     interaction test ([siga_screen()]).
#'   \item **Genetic algorithm search** ([siga_search()]). Every individual of
#'     the population is a vector of real-valued genes, one per candidate
#'     biomarker, that encodes whether the biomarker is used and which threshold
#'     or set of levels defines the subgroup. Individuals are scored by the
#'     p-value of the treatment-by-subgroup interaction and evolve with the
#'     operators of the \pkg{gatbxr} toolbox (ranking, stochastic universal
#'     sampling, line recombination and breeder genetic algorithm mutation) in
#'     several isolated sub-populations.
#'   \item **Resampling**: the search is repeated on random 63.2% subsamples
#'     and the subgroups found are evaluated on the held-out 36.8%. The average
#'     held-out treatment effect controls multiplicity (parameter `eta`) and the
#'     frequency with which every biomarker is selected controls complexity
#'     (parameter `theta`).
#'   \item **Final search** on the complete data with the frequently selected
#'     biomarkers.
#' }
#' The complete procedure is run by [siga()].
#'
#' @section Main functions:
#' \describe{
#'   \item{[siga()]}{the complete SIGA procedure.}
#'   \item{[update.siga()]}{re-apply the decision rule with other values of
#'     `eta` and `theta` without repeating the resampling.}
#'   \item{[siga_search()]}{the (multi-population) genetic algorithm search.}
#'   \item{[siga_screen()]}{candidate biomarker prescreening.}
#'   \item{[siga_cv()]}{cross-validation of the complete procedure.}
#'   \item{[subgroup_effect()]}{treatment effects in a subgroup, its complement
#'     and the treatment-by-subgroup interaction.}
#'   \item{[siga_simulate()]}{data from the simulation models of the SIGA paper.}
#' }
#'
#' @references
#' Zhao L, Zhang W, Wu Y, Cao L, Wang L, Li K (2024). A promising subgroup
#' identification method based on a genetic algorithm for censored survival
#' data. *Journal of Biopharmaceutical Statistics* 34(1), 55--77.
#' \doi{10.1080/10543406.2023.2170397}
#'
#' Chipperfield A, Fleming P, Pohlheim H, Fonseca C (1994). Genetic Algorithm
#' Toolbox for use with MATLAB. University of Sheffield.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom survival Surv
#' @export
survival::Surv

## The gatbxr operators are called by name; gatbxr::recombin() looks the name
## up from the calling environment, so the operators are imported here.
#' @importFrom gatbxr mutbga recdis recint reclin rws sus
## usethis namespace: start
#' @importFrom stats median pchisq pf pnorm predict pt qlogis qnorm quantile
#'   rbinom rexp rnorm runif setNames terms update
#' @importFrom graphics abline axis barplot box hist legend lines par plot
#'   points segments text title
#' @importFrom grDevices dev.flush dev.hold gray
#' @importFrom utils combn head
## usethis namespace: end
NULL
