#' @exportS3Method
print.downscalr = function(x,lu.from=NULL, lu.to=NULL, times=NULL, value=NULL, downscale.value=NULL, target=NULL, ...){

  check.targets =
    x$ds.inputs$targets %>%
      dplyr::group_by(lu.from,lu.to,times) %>%
      dplyr::summarise(target=sum(value), .groups = "keep") %>%
      dplyr::left_join(
            x$out.res %>%
              dplyr::group_by(lu.from,lu.to,times) %>%
              dplyr::summarise(downscale.value=sum(value), .groups = "keep"),
            by = c("lu.from", "lu.to", "times")) %>%
      dplyr::mutate(differece=(downscale.value-target)) %>%
      dplyr::ungroup()

  print(check.targets)
}

