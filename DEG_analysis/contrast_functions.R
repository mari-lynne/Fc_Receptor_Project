#Just contrast matrix function ####
make_contrasts <- function (group, control, delim="_vs_", des_mat){
  #/ define groups and baseline to make contrasts
  
  suppressMessages(require(limma))
  
  #Checks
  if(is.null(group)) stop("Error: group arg is missing")
  
  #/ ensure unique group levels
  group <- sort(unique(as.character(group)))
  
  #Write limma code by pasting groups
  #/ if control var is present, compare all groups to control
  #/  else make all comparisons using combn function
  if (!missing(control)){
    combo <- paste0(group,"-",  control)
  } else{
    combo <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])})
  }
  
  #/ make contrasts
  if (!missing(des_mat)){
  contrasts<- limma::makeContrasts(contrasts=combo, levels=colnames(des_mat))
  }else{
    contrasts<- limma::makeContrasts(contrasts=combo, levels=group)
    message("No Design Matrix provided, using only defined contrasts for matrix")
  }
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  #Todo: levels need to be design matrix, poss modify to have group or design option later
  #Only do lmfit step if design and fit are supplied
  
  return(contrasts)
}

#output testing
#results <- make_contrasts(group = c("V7", "V1", "D0", "V0"), control = c("V0"), des_mat = design)

#Combo with ebayes #####
contrast_2_lm <- function(group, control, delim="_vs_", des_mat, efit, topTab="TRUE"){
  #Define groups and baseline to make contrasts
  #Input design and previous model fit model for it to work
  #Toptab gives option to output results
  
  #/ Checks
  suppressMessages(require(limma))
  if(is.null(group)) stop("Error: group arg is missing")
  
  #/ Ensure unique group levels
  group <- sort(unique(as.character(group)))
  
  #/ Define contrasts 
  if (!missing(control)){ #compare to control
    combo <- paste0(group,"-",  control)
  } else{ #make all pairwise comparisons
    combo <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])}) 
  }
  
  #/ Make contrast matrix
  if (!missing(des_mat)){
    contrasts<- limma::makeContrasts(contrasts=combo, levels=colnames(design))
  }else{ #No design + no efit
    contrasts<- limma::makeContrasts(contrasts=combo, levels=group)
    message("No Design Matrix provided, using only defined contrasts for matrix")
  } 
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  
 
  #/ Model fit and deg results
  if (!missing(efit) & !missing(des_mat) & topTab == "TRUE"){ #Requires efit and dmat args
    fit2 <- contrasts.fit(efit, contrasts)
    fit2 <- eBayes(fit2)
    message("Performing ebayes fit of linear model")
    top_results <- list() #TopTable Results
    for (i in colnames(contrasts)){
      top_results[[i]] <- topTable(fit2, coef = i, number = Inf)
      limma_list <- list(contrasts = contrasts, fit2 = fit2, top_results = top_results)
    }
    message("Contrast matrix (contrasts), ebayes (fit2), top DEGs (top_results) saved in list\n
          Subset list with either $ or [[]] for results")
  } else if (!missing(efit) & !missing(des_mat) & topTab != "TRUE"){ #Without toptable option
    fit2 <- contrasts.fit(efit, contrasts)
    fit2 <- eBayes(fit2)
    message("Performing ebayes fit of linear model")
    limma_list <- list(contrasts = contrasts, fit2 = fit2)
  } else if (missing(efit) | missing(des_mat)) {
    limma_list <- list(contrasts)
    warning("No linear model or design matrix supplied, returning contrast matrix only")
  }
  return(limma_list)
}

#Output testing 
results2 <- contrast_2_lm(group = c("V7", "V1", "D0"), control = c("V0"), efit=fit, des_mat = design)

