#' @exportS3Method summary downscalr
summary.downscalr = function(object, lu.from=NULL, lu.to=NULL, times=NULL, value=NULL, downscale.value=NULL, target=NULL, difference= NULL, miss.threshold=0.01, ...){

  summarylist <- list()

  check.targets =
    object$ds.inputs$targets %>%
    dplyr::group_by(lu.from,lu.to,times) %>%
    dplyr::summarise(target=sum(value), .groups = "keep") %>%
    dplyr::left_join(
      object$out.res %>%
        dplyr::group_by(lu.from,lu.to,times) %>%
        dplyr::summarise(downscale.value=sum(value), .groups = "keep"),
      by = c("lu.from", "lu.to", "times")) %>%
    dplyr::mutate(difference=ifelse(is.na((downscale.value-target)/target),0,(downscale.value-target)/target)) %>%
    dplyr::ungroup()

  summarylist[['check.targets']] <- check.targets

  summarylist[['max.miss']] <- check.targets %>% filter(abs(difference)==max(abs(difference)))

  summarylist[['table.miss']] <- check.targets %>% filter(abs(difference)>=miss.threshold)

  summarylist[['miss.threshold']] <- miss.threshold

  return(summarylist)

}


#' @exportS3Method print downscalr
print.downscalr = function(x, ...){
  to.print <- summary.downscalr(x)

  cat('Largest deviation from target in per cent of target:')
  cat('\n')
  print(to.print[["max.miss"]])
  cat('\n')
  if(nrow(to.print[["table.miss"]])!=0){
    cat(paste0('Transitions with deviation from target larger or equal to ',to.print[['miss.threshold']],' per cent:'))
    cat('\n')
    print(to.print[["table.miss"]])
  } else {
    cat(paste0('No deviation from target larger than or equal to ',to.print[['miss.threshold']],' per cent:'))
  }

}


# to add is a level map in the summary object and print aggregates in print output 123
