# =============================================================================
# Author: Haodong Tian
# Description: Helper function to convert GenomicSEM output objects into
#   semPlot-compatible model objects for path diagram visualization.
# Note: Comments and sensitive paths have been cleaned for journal code review.
# =============================================================================

## Reference: https://groups.google.com/g/genomic-sem-users/c/h250OMkp_CM/m/S07GtuPYAAAJ

create_semPlotModel <- function(gsem.object=GWISoutput, est.label="STD_All"){

  object <- gsem.object$results
  object$free <- 0
  numb <- 1:length(which(object$op != "~~"))
  object$free[which(object$op != "~~")] <- numb
  varNames  <- lavaanNames(object, type = "ov")
  factNames <- lavaanNames(object, type = "lv")
  factNames <- factNames[!factNames %in% varNames]
  n <- length(varNames)
  k <- length(factNames)
  if (is.null(object$label))
    object$label <- rep("", nrow(object))
  semModel <- new("semPlotModel")
  object$est <- object[, est.label]
  if (is.null(object$group))
    object$group <- ""
  semModel@Pars <- data.frame(label = object$label,
                              lhs   = ifelse(object$op == "~" | object$op == "~1", object$rhs, object$lhs),
                              edge  = "--",
                              rhs   = ifelse(object$op == "~" | object$op == "~1", object$lhs, object$rhs),
                              est   = object$est,
                              std   = NA,
                              group = object$group,
                              fixed = object$free == 0,
                              par   = object$free,
                              stringsAsFactors = FALSE)
  semModel@Pars$edge[object$op == "~~"]               <- "<->"
  semModel@Pars$edge[object$op == "~*~"]              <- "<->"
  semModel@Pars$edge[object$op == "~"]                <- "~>"
  semModel@Pars$edge[object$op == "=~"]               <- "->"
  semModel@Pars$edge[object$op == "~1"]               <- "int"
  semModel@Pars$edge[grepl("\\|", object$op)]         <- "|"
  semModel@Thresholds <- semModel@Pars[grepl("\\|", semModel@Pars$edge), -(3:4)]
  semModel@Pars       <- semModel@Pars[!object$op %in% c(":=", "<", ">", "==", "|", "<", ">"), ]
  semModel@Vars <- data.frame(name     = c(varNames, factNames),
                              manifest = c(varNames, factNames) %in% varNames,
                              exogenous = NA,
                              stringsAsFactors = FALSE)
  semModel@ObsCovs  <- list()
  semModel@ImpCovs  <- list()
  semModel@Computed <- FALSE
  semModel@Original <- list(object)
  return(semModel)
}


# Example usage:
# semPaths(create_semPlotModel(GWISoutput), layout="tree2", whatLabels="est")
