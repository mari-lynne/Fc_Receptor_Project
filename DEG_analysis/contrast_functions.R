#Just cm function ####

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
results <- make_contrasts(group = c("V7", "V1", "D0", "V0"), control = c("V0"), des_mat = design)



#Combo with ebayes #####


contrast_2_lm <- function(group, control, delim="_vs_", des_mat, efit, topTab="TRUE"){
  #/ define groups and baseline to make contrasts
  #/ input design and previous model fit model for it to work
  #Toptab gives option to output results
  
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
    contrasts<- limma::makeContrasts(contrasts=combo, levels=colnames(design))
  }else{
    contrasts<- limma::makeContrasts(contrasts=combo, levels=group)
    message("No Design Matrix provided, using only defined contrasts for matrix")
  }
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
 
  #/model fit
  if (!missing(efit) & !missing(des_mat)){ 
    fit2 <- contrasts.fit(efit, contrasts)
    fit2 <- eBayes(fit2)
    message("Performing ebayes fit of linear model")
    limma_list <- list(contrasts = contrasts, fit2 = fit2)
  } else { #Only run ebayes model if initial lm and design are supplied
    limma_list <- list(contrasts)
    warning("No linear model supplied, returning contrast matrix only")
  }
  
  #/Toptable
  if (topTab == "TRUE" & !missing(efit) & !missing(des_mat)){
    top_results <- list()
    for (i in colnames(cm)){
      top_results[[i]] <- topTable(fit2, coef = i, number = Inf)
      limma_list <- list(contrasts = contrasts, fit2 = fit2, top_results = top_results)
    }
    message("Contrast matrix (contrasts), ebayes (fit2), top DEGs (top_results) saved in list\n
          Subset list with either $ or [[]] for results")
  } else { #No top table results #else if maybe necessary
    limma_list <- list(contrasts, efit)
    message("Contrast matrix (contrasts), ebayes (fit2) saved in list.\n
          Subset list with either $ or [[]] for results")
  }
  
  return(limma_list)
}
  
#Output testing ####
results2 <- contrast_2_lm(group = c("V7", "V1", "D0"), control = c("V0"),efit=fit, des_mat = design)

