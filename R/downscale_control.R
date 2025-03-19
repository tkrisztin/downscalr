#' Set Options for Downscaling Solver
#'
#' This function configures the downscalR solver with customizable options. It allows the selection of the solver function and adjustment of solver behavior and stopping conditions, tailored for specific downscaling needs.
#'
#' @param solve_fun The downscaling function to use. Valid options are "\code{solve_biascorr}" for bias correction and "\code{solve_notarget}" for non-targeted downscaling. Defaults to "\code{solve_biascorr}".
#' @param algorithm Specifies the solver algorithm. The default algorithm is "\code{NLOPT_LN_SBPLX}". For more information on available algorithms, refer to the \code{\link[nloptr]{nloptr}} package.
#' @param xtol_rel The relative tolerance level of the solver, controlling solution precision in relative terms. Default is \code{1.0e-20}, where a lower value indicates higher precision.
#' @param xtol_abs The absolute tolerance level of the solver, specifying solution precision in absolute terms. Default is \code{1.0e-20}, with a lower value indicating higher precision.
#' @param maxeval The maximum number of evaluations the solver will perform, limiting computational effort. Default is \code{1600}.
#' @param MAX_EXP Sets the numerical cutoff to prevent extremely high values in the Multinomial Logit (MNL) function. Default is \code{log(.Machine$double.xmax)}.
#' @param cutoff An optional cutoff in the share of a pixel to avoid MNL values close to zero. Default is \code{0}.
#' @param max_diff The threshold for acceptable difference to targets. If exceeded, the estimation is redone. Default is \code{1.0e-8}.
#' @param ref_class_adjust_threshold The threshold for adjusting the reference class in the MNL. If the implied target share for the reference class is below this, additional areas are created for numerical stability and removed after solving. Default is \code{1.0e-8}.
#' @param err.txt An optional error message for debugging purposes. Empty by default.
#'
#' @return Returns a list with the specified options for the downscaling solver.
#'
#' @details
#' The function returns a list of parameters to be used by the solver. It is essential to provide compatible parameters to ensure accurate downscaling results. The function will stop with an error message if \code{solve_fun} is not correctly specified.
#'
#' @export downscale_control
#'
#' @examples
#' # Default settings
#' opts_default = downscale_control()
#'
#' # Custom settings
#' opts_custom = downscale_control(algorithm = "some_other_algorithm",
#'                                 xtol_rel = 1e-5, xtol_abs = 1e-5,
#'                                 maxeval = 2000)
downscale_control = function(solve_fun = "solve_biascorr",algorithm = "NLOPT_LN_SBPLX",
                                  xtol_rel = 1.0e-20,xtol_abs = 1.0e-20,maxeval = 1600,
                                  MAX_EXP = log(.Machine$double.xmax),cutoff = 0,
                                  max_diff = 1.0e-8,ref_class_adjust_threshold = 1.0e-8,
                             err.txt = "") {
  if (!solve_fun %in% c("solve_biascorr","solve_notarget")) {stop("solve_fun not correctly specified.")}
  return(list(solve_fun = solve_fun,algorithm = algorithm,xtol_rel = xtol_rel,xtol_abs = xtol_abs,
              maxeval = maxeval,MAX_EXP = MAX_EXP,cutoff = cutoff,
              max_diff = max_diff,ref_class_adjust_threshold=ref_class_adjust_threshold,
              err.txt = err.txt
  ))
}


#' Set options for population downscaling
#'
#' This function can be called to specify additional options to the downscalR solver.
#' The user can change the solver function, specify solver behavior, and stopping conditions.
#'
#' @param solve_fun The name of the downscaling function to use. Has to be a
#' valid R function. Defaults to \code{"solve_biascorr"}.
#' @param MAX_EXP Numerical cutoff for MNL function.
#' @param err.txt Error text for caller identification (used for debugging)
#' @param max_diff If difference to targets is larger, redo the estimation (helps to avoid convergence errors)
#'
#' @return List with default options for bias correction solver
#'
#' @export downscale_control_pop
#'
#' @examples
#' opts1 = downscale_control_pop()
downscale_control_pop = function(solve_fun = "solve_biascorr",
                                 MAX_EXP = log(.Machine$double.xmax),
                                 max_diff = 1.0e-8,
                             err.txt = "") {
  return(list(solve_fun = solve_fun,
              MAX_EXP = MAX_EXP,
              max_diff = max_diff,
              err.txt = err.txt
  ))
}
