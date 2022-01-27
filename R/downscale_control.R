#' Get default options for bias corrections solver
#'
#' @param algorithm Solver algorithm (see the \code{\link[nloptr]{nloptr}} package for documentation and more detail).
#' @param xtol_rel Relative tolerance of solver.
#' @param xtol_abs Absolute solver tolerance.
#' @param maxeval Maximum evaluation of solver.
#' @param MAX_EXP Numerical cutoff for MNL function.
#' @param cutoff Optional cutoff to avoid MNL values close to zero.
#' @param err.txt Error text for caller identification (used for debugging)
#' @param max_diff If difference to targets is larger, redo the estimation (helps to avoid convergence errors)
#' @param redo Maximum number of repeats
#'
#' @return List with default options for bias correction solver
#'
#' @details Call this function if you want to change default options for the bias corrections solver.
#'
#' @export downscale_control
#'
#' @examples
#' opts1 = downscale_control()
downscale_control = function(algorithm = "NLOPT_LN_SBPLX",
                                  xtol_rel = 1.0e-20,xtol_abs = 1.0e-20,maxeval = 1600,
                                  MAX_EXP = log(.Machine$double.xmax),cutoff = 0,
                                  redo = 2,max_diff = 10^-8,err.txt = "") {
  return(list(algorithm = algorithm,xtol_rel = xtol_rel,xtol_abs = xtol_abs,
              maxeval = maxeval,MAX_EXP = MAX_EXP,cutoff = cutoff,redo = redo,
              max_diff = max_diff,err.txt = err.txt
  ))
}
