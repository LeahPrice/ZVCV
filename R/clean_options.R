# Add defaults where option fields are missing and check reasonable polynomial order
clean_options <- function(options, N, d_full){
  
  if (any(c("polyorder","regul_reg","alpha_elnet","nfolds","apriori","intercept","polyorder_max") %in% names(options))){
    options <- rep(list(options),1)
  }
  num_options <- length(options)
  
  for (i in 1:num_options){
    defined_terms <- length(options[[i]])
    correctly_defined_terms <- 0
    if (!("polyorder" %in% names(options[[i]]))) { options[[i]]$polyorder <- 2 } else { correctly_defined_terms <- correctly_defined_terms + 1 }
    if (!("regul_reg" %in% names(options[[i]]))) { options[[i]]$regul_reg <- TRUE } else { correctly_defined_terms <- correctly_defined_terms + 1 }
    if (!("alpha_elnet" %in% names(options[[i]]))) { options[[i]]$alpha_elnet <- 1 } else { correctly_defined_terms <- correctly_defined_terms + 1 }
    if (!("nfolds" %in% names(options[[i]]))) { options[[i]]$nfolds <- 10 } else { correctly_defined_terms <- correctly_defined_terms + 1 }
    if (!("apriori" %in% names(options[[i]]))) { options[[i]]$apriori <- 1:d_full } else { correctly_defined_terms <- correctly_defined_terms + 1 }
    if (!("intercept" %in% names(options[[i]]))) { options[[i]]$intercept <- TRUE} else { correctly_defined_terms <- correctly_defined_terms + 1 }
    if ("polyorder_max" %in% names(options[[i]])){ correctly_defined_terms <- correctly_defined_terms + 1 }
    
    if (correctly_defined_terms != defined_terms){
      stop("Check options. It should only include the terms polyorder, regul_reg, alpha_elnet, nfolds, apriori, intercept or polyorder_max")
    }
    
    d <- length(options[[i]]$apriori)
    if (d==1 & !("polyorder_max" %in% names(options[[i]]))){
      options[[i]]$polyorder_max <- Inf
    } else if (d > 1){ # Check for potentially large design matrices when d>1
      temp_max_polyorder <- 1 # identify the largest polynomial order such that the regression design matrix has no more than 10^6 elements
      while (N*choose(d + temp_max_polyorder,d)<10^7){ 
        temp_max_polyorder <- temp_max_polyorder + 1
      }
      
      if (is.infinite(options[[i]]$polyorder) & !("polyorder_max" %in% names(options[[i]]))) { # Tell users about maximum polynomial orders
        options[[i]]$polyorder_max <- temp_max_polyorder
        warning(sprintf("To prevent memory issues, the maximum polynomial order for option %d has been set to %d.\n",i,temp_max_polyorder))
      } else if (!is.infinite(options[[i]]$polyorder) & temp_max_polyorder < options[[i]]$polyorder){
        input_continue <- readline(prompt=sprintf("\nThe polynomial order for option %d will result in a design matrix\n of size %d by %d. This may require long computation time\n or large amounts of memory. Are you sure you wish to\n continue? Y/N: ",i,N,choose(d+options[[i]]$polyorder,d)))
        input_continue <- as.character(input_continue)
        if (toupper(input_continue) == "N"){
          stop("To fix this issue, consider reducing the polynomial order or using a subset through the apriori option.")
        } else if (toupper(input_continue) != "Y"){
          stop("Input must be either Y or N.")
        }                
      }
    }
  }
  return(options)
}

