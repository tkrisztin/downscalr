#' Iterated Grid Search for Function Optimization with Additional Parameters
#'
#' Performs an iterated grid search over a given range to find the parameter that minimizes a provided function.
#' The search starts by evaluating the function over a sequence of points between a specified minimum and maximum value.
#' The search range is iteratively refined around the minimum point until the maximum number of iterations is reached
#' or the desired precision is achieved. Optionally, the search grid can be transformed exponentially.
#'
#' @param min_param Numeric. The minimum value of the parameter range to start the search.
#' @param max_param Numeric. The maximum value of the parameter range to start the search.
#' @param func Function. The function to be minimized. It should accept a numeric parameter and return a numeric value.
#' @param step_length Integer. The number of points to evaluate in each iteration. Default is 100.
#' @param max_iterations Integer. The maximum number of iterations to perform. Default is 5.
#' @param precision_threshold Numeric. The precision threshold for stopping the search. Once the difference between the
#'                             min and max of the refined range is less than this value, the search will stop. Default is 1e-3.
#' @param exp_transform Logical. If TRUE, the search grid will be transformed exponentially (on the log scale). Default is FALSE.
#' @param ... Additional parameters to be passed to the function `func`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{best_param}{The parameter that gives the minimum function value.}
#'   \item{best_value}{The minimum value of the function.}
#' }
#'
#' @examples
#' # Define a sample function to minimize
#' test_function <- function(x, a = 1) { a * (x - 3)^2 + 2 }
#'
#' # Perform the iterated grid search with additional parameter and exponential transform
#' result <- iterated_grid_search(min_param = 0, max_param = 10, func = test_function, a = 2, exp_transform = TRUE)
#' print(result)
#'
#' @export
iterated_grid_search <- function(min_param,
                                 max_param,
                                 func,
                                 step_length = 100,
                                 max_iterations = 5,
                                 precision_threshold = 1e-3,
                                 exp_transform = FALSE,
                                 ...) {
  for (iteration in 1:max_iterations) {
    # Create a sequence between min_param and max_param with length equal to step_length
    search_grid <- seq(min_param, max_param, length.out = step_length)

    # Optionally transform the grid exponentially
    if (exp_transform) {
      search_grid <- exp(search_grid)
    }

    # Apply the function to each value in the search_grid, passing on additional arguments
    func_values <- sapply(search_grid, func, ...)

    # Find the index of the minimum value of the function
    min_index <- which.min(func_values)

    # Get the corresponding minimum value from the grid
    best_param <- search_grid[min_index]
    best_value <- func_values[min_index]

    # Calculate the new range: left and right of the minimum
    if (min_index == 1) {
      new_min_param <- search_grid[min_index]
      new_max_param <- search_grid[min_index + 1]
    } else if (min_index == length(search_grid)) {
      new_min_param <- search_grid[min_index - 1]
      new_max_param <- search_grid[min_index]
    } else {
      new_min_param <- search_grid[min_index - 1]
      new_max_param <- search_grid[min_index + 1]
    }

    # Check if the range is within the precision threshold
    if (abs(new_max_param - new_min_param) < precision_threshold) {
      return(list(best_param = best_param, best_value = best_value))
    }

    # Update the min and max parameters for the next iteration
    min_param <- ifelse(exp_transform,log(new_min_param),new_min_param)
    max_param <- ifelse(exp_transform,log(new_max_param),new_max_param)
  }

  # Return the best parameter and the corresponding function value after max iterations
  return(list(best_param = best_param, best_value = best_value))
}
