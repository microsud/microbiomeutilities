
# Check group variable
.check.group <- function(x, group) {
  
  if(isFALSE(any(group %in% sample_variables(x)))){
    stop("'group' variable must be in sample_variables(x)")
  } 
  
}

