
##
##
##
##
########################################################################################################################################################
########################################################################################################################################################
##
## Compute cohort-specific absolute impacts
##
########################################################################################################################################################
########################################################################################################################################################
##
##
##
##
##
##
##

#' Network propagation module: Compute absolute cohort-specific impact matrix
#'
#' This function computes an absolute cohort-specific impact matrix for a given data set (data) under consideration of a given network (networkName) using network propagation. Only significant regulatory links (pValCutoff) that are not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene are included in the computations. First, a basic initial cohort-specific impact matrix is computed. Next, this matrix is used to iteratively compute the absolute cohort-specific impact matrix stopping if the sum of the column sums of the current impact matrix increases less than a given threshold (colSumsThreshold) in comparison to the previous impact matrix. The basic inital impact matrix and the resulting final impact matrix are saved using pre-defined file naming conventions. The resulting final impact matrix can be used to determine absolute impacts of source on target genes. Absolute impacts allow to identify those genes that have the greatest impact on other genes. Absolute impacts do not distinguish between potential activator and inhibitory contributions along the network paths.
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation of impact matrix if improvement is less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlow_CohortSpecificAbsoluteImpacts}}, \code{\link{loadNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}, \code{\link{getAverageImpacts_CohortSpecificAbsoluteImpacts}}, \code{\link{getImpacts_CohortSpecificAbsoluteImpacts}}, \code{\link{getRawImpacts_CohortSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}}, ...
#' @return Saves basic initial and resulting final impact matrices. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts( data = data, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath )
#' 
computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts <- function( data, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
    ##Compute basic network flow matrix
    computeBasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts( data = data, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, path = path, output = output )
    
    ##Compute final network flow matrix
    computeNetworkFlow_CohortSpecificAbsoluteImpacts( dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output, patientSpecific = FALSE )    
}



#' Network propagation module: Compute basic initial absolute cohort-specific impact matrix
#'
#' This function computes a basic initial cohort-specific impact matrix that is required for the iterative computation of an absolute cohort-specific impact matrix. The basic matrix is computed for a given data set (data) under consideration of a given network (networkName) using network propagation. Only significant regulatory links (pValCutoff) that are not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene are included in the computations. The basic inital impact matrix is saved using pre-defined file naming conventions. This matrix can be used to compute an absolute cohort-specific impact matrix using the function \code{\link{computeNetworkFlow_CohortSpecificAbsoluteImpacts}}. This function allows to separate the computation of the initial from the final impact matrix.
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlow_CohortSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}.
#' @return Saves basic initial impact matrix. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' computeBasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts( data = data, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#' 
computeBasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts <- function( data, dataSetName, networkName, pValCutoff, localGeneCutoff, path, output = TRUE )
{
  ##Load network and remove non-significant predictors
  network = loadNetworkWithFilteringForSignificantPredictors( networkName = networkName, pValCutoff = pValCutoff, path = path, output = output )

  ##Check if genes in the data set considered for prediction are in the same order of genes as in the data set used to learn the network
  if( !identical( data$genes, network$responseGenes ) )
  {
      print( paste( "Error: Genes in the given data set '", dataSetName, "' are not in the same order as the genes in the data set used to learn the network '", networkName, "'.", sep = "" ) )
      print( "" )
      print( " - Please ensure that the given data set has exactly the same order of genes." )
      print( " - Remove genes from the given data set if they were not contained in the data set used for network learning." )
      print( " - Add missing genes with artifical measurement values of zero to the given data set." )
      
      return()
  }
  
  ##
  ##Basic variables
  ##    
  #Number of response variables: all genes in the data set
  nrResponses = length( data$genes )

  #Number of predictor variables (gene-specific copy number + expression levels)
  nrPredictors = nrow( data$Y ) + 1

  #Number of experiments
  T = ncol( data$U )  
  
  ##
  ##Basic matrices: Y and U contain data of considered response genes; M must be complete
  ##   
  #Response data
  #Gene expression levels of response genes: Y matrix -> number of experiments x number of response genes
  Y   = t( data.matrix( data$Y ) )
  
  #Predictor data
  #Copy numbers of response genes: U matrix -> number of experiments x number of response genes
  U   = t( data.matrix( data$U ) )      
  #Gene expression levels of predictor genes: M matrix -> number of experiments x number of predictor genes
  M   = t( data.matrix( data$Y ) )
  
  
 
  ##Network flow matrix
  #- rows: predictor genes
  #- colums: response genes
  #- rows and columns are exchanged compared to the network matrix to enable useful matrix multiplications
  #- square matrix: nrResponses x nrResponses to enable network flow analyses via matrix multiplication
  #
  
  networkFlowMatrix = matrix( 0, nrow = nrResponses, ncol = nrResponses )
  colnames( networkFlowMatrix ) <- data$genes
  rownames( networkFlowMatrix ) <- data$genes

  ##Direct copy number contributions
  directCopyNumberContribution = matrix( 0, nrow = nrResponses, ncol = 1 )
  rownames( directCopyNumberContribution ) <- data$genes
  colnames( directCopyNumberContribution ) <- "CopyNumberContribution"
  
  ##Explained variance
  correlationStatistics = matrix( 0, nrow = nrResponses, ncol = 1 )
  rownames( correlationStatistics ) <- data$genes
  colnames( correlationStatistics ) <- "AverageCorrelation"
  

  ##Explained variance
  explainedVariance = matrix( 0, nrow = nrResponses, ncol = 1 )
  rownames( explainedVariance ) <- data$genes
  colnames( explainedVariance ) <- "ExplainedVariance"

  if( output )
  {
      print( "Compute cohort-specific basic network flow matrix" )
  }
  
  for( i in 1:nrResponses )
  {
      ##
      ##Response vector
      ##
      y = Y[ , i ]
      
      ##
      ##Remove predictors that are on the same chromosome like respone gene i in distance <= localGeneCutoff
      ##
      Mdummy = M
      
      if( localGeneCutoff > 0 )
      {
        currentResponseGene    = data$genes[ i ]
        currentResponseGeneChr = data$chr[ i ]

        #Determine downstream predictors
        removePredictorPos = c()
        for( v in max( 1, i - localGeneCutoff ):max( 1, i - 1 ) )
        {
          localNeighborGeneChr = data$chr[ v ]
  
          #If local neighbor is still on same chromosome then mark it for later removement
          if( localNeighborGeneChr == currentResponseGeneChr )
          {
            removePredictorPos = c( removePredictorPos, v )
          }      
        }

        #Determine upstream predictors
        for( v in min( nrResponses, i + 1 ):min( nrResponses, i + localGeneCutoff ) )
        {
          localNeighborGeneChr = data$chr[ v ]
    
          #If local neighbor is still on same chromosome then mark it for later removement
          if( localNeighborGeneChr == currentResponseGeneChr )
          {
            removePredictorPos = c( removePredictorPos, v )
          }      
        }

        #Remove predictors (note that the response gene itself is already removed during the network computation -> all network weights are zero for the response gene)
        Mdummy[ , removePredictorPos ] = 0      
      }

      ##
      ##Predictor matrix
      ##
      #Add copy number of response gene as a predictor
      specificM = cbind( U[ , i ], Mdummy )
                
      
      ##
      ##Compute contributions of individual predictors on the expression level of response gene i
      ##      
      #Contributions of individual predictors (first entry is the copy number effect the remaining entries reflect effects from trans predictors)
      predictorContributions = computePredictorProportions_CohortSpecificAbsoluteImpacts( data = specificM, regressionParameters = network$G[ i,  ], predictorGenes = network$predictors, responseGenes = network$responseGenes )

    
      ##
      ##Compute the correlation between predicted and observed gene expression levels of gene i
      ##
      #Predict expression level ignoring offset (later correlation evaluation is invariant to the offset)
      Yhat = specificM %*%  network$G[ i,  ]
      
      corStatistics = 0
      if( var( Yhat, na.rm = TRUE ) == 0 )
      {
        if( output )
        {
          print( paste( data$genes[ i ], "not predictable" ) )
        }
      }
      else
      {
        #Compute correlation and p-value
        dummy = cor.test( y, Yhat, alternative = "greater" )

        corStatistics = dummy[[ 4 ]] #Correlation
      }
      
      ##
      ##Compute the proportion of variance for gene i that is explained by the network
      ##            
      #Prediction not possible: e.g. no predictors available
      if( is.na( corStatistics ) )
      {
        corStatistics = 0
      }
      
      #Network cannot predict the expression level correctly
      if( corStatistics < 0 )
      {
        corStatistics = 0 
      }
      
      #Proportion of variance of gene i explained by the network
      squareCorrelation = corStatistics^2
      
      ##Multiply the predictor contributions with the explained variance
      predictorContributions = predictorContributions * squareCorrelation
      
      ##Transfer predictor contributions to the corresponding network flow column
      networkFlowMatrix[ , i ] = predictorContributions[ 2:( nrResponses + 1 ) ]
      
      ##
      ##Integrate direct copy number effect
      ##
      #- Entry [i,i] is always zero, since no gene can have a trans-contribution on its own expression
      #- We can use [i,i] to integrate the gene-specific copy number contribution 
      networkFlowMatrix[ i, i ] = predictorContributions[ 1 ]

      ##Transfer direct copy number contribution
      directCopyNumberContribution[ i, 1 ] = predictorContributions[ 1 ]
      
      ##Transfer average correlation
      correlationStatistics[ i, 1 ] = corStatistics
      
      ##Transfer explained variance
      explainedVariance[ i, 1 ] = squareCorrelation
  }
  

  ##
  ##Save basic network flow matrix
  ##
  if( output )
  {
      print( "Save output:" )
  }

  subPath = "/NetworkPropagation/NetworkFlow/BasicNetworkFlowMatrices/"
  outputFile =  paste( path, subPath, dataSetName, "_BasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".Rout", sep = "" )
  
  dummy <- Matrix( networkFlowMatrix, sparse = TRUE )
  networkFlowMatrix <- dummy
  
  save( networkFlowMatrix, directCopyNumberContribution, explainedVariance, correlationStatistics, file = outputFile )
  
  if( output )
  {
      print( outputFile )
  }

}



#
#Compute the proportion that individual predictors have on the prediction of the expression level of a gene: internal function.
#
computePredictorProportions_CohortSpecificAbsoluteImpacts <- function( data, regressionParameters, predictorGenes, responseGenes )
{
    activePredictors = which( regressionParameters != 0 )    
    N = length( activePredictors )
    
    predictorProportions = rep( 0, N )
    
    #Global return values
    copyNumberContribution = 0 
    globalTransPredictorContributionVector = rep( 0, times = length( responseGenes ) )
    
    ##
    ##Compute absolute proportion of individual predictors
    ##
    if( N > 0 )
    {
      ##Compute proportion
      for( i in 1:N )
      {
        #Average absolute contribution of the predictor is obtained by
        #- multiplying the measurements of the predictor with its corresponding regression parameter
        #- followed by averaging over the different tumor samples
        predictorProportions[ i ] = mean( abs( data[ , activePredictors[ i ] ] * regressionParameters[ activePredictors[ i ] ] ), na.rm = TRUE )
      }
      
      #Standardization
      predictorProportions = predictorProportions / sum( predictorProportions )
      
      ##
      ##Create global predictor contribution vector
      ##
      #Determine copy number contribution if active and remove copy number predictor from the list of active predictors, because we are later interested in the influence of trans effects         
      if( activePredictors[ 1 ] == 1 )
      {
        copyNumberContribution = predictorProportions[ 1 ]  
        activePredictors = setdiff( activePredictors, 1 )
        if( N > 1 )
        {
          #remove copy number proportion
          predictorProportions = predictorProportions[ 2:N ]
        }
      }
      
      ##This shift is necessary, since the regression parameters start with the copy number predictor
      activePredictors = activePredictors - 1      
      
      N = length( activePredictors ) #if there is only a copy number predictor, N will be zero
      if( N > 0 )
      {
        activePredictorGenes = responseGenes[ activePredictors ]
  
        mapPos = match( activePredictorGenes, responseGenes )
  
        ##Map predictor genes to their corresponding position in the trans-contribution vector (order: responseGenes)
        dummy = which( !is.na( mapPos ) )
        uniqueMapPos  = mapPos[ dummy ]  
        globalTransPredictorContributionVector[ uniqueMapPos ] = predictorProportions[ dummy ]
      }      
      
    }
    
    res  = checkContributions( res = c( copyNumberContribution, globalTransPredictorContributionVector ) )
   
    #First entry is the copy number contribution the remaining entries reflect the trans predictors
    return( res )
} 



#
#Remove NaN, NAs and Inf if present. Such entries are replaced by zero to guarantee a proper network flow computation: internal function.
#
#- NaN, NA, Inf can occur if no predictors of a response gene have been measured (e.g. genes not present in a test data set and one has added the predictor genes to the test data set to be able to apply the learned network)
#
checkContributions <- function( res )
{
    cand = which( is.nan( res ) )
    if( length( cand ) > 0 )
    {
      res[ cand ] = 0
    }
    
    cand = which( is.na( res ) )
    if( length( cand ) > 0 )
    {
      res[ cand ] = 0
    }
    
    cand = which( is.infinite( res ) )
    if( length( cand ) > 0 )
    {
      res[ cand ] = 0
    }
    
    #First entry is the copy number contribution the remaining entries reflect the trans predictors
    return( res )
}



#' Network propagation module: Compute absolute cohort-specific impact matrix
#'
#' This function computes an absolute cohort-specific impact matrix for a given data set (dataSetName) under consideration of a given network (networkName) using network propagation. The basis of this computation is the pre-computed basic initial cohort-specific impact matrix that only considered significant regulatory links (pValCutoff) that were not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene (\code{\link{computeBasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}). The absolute cohort-specific impact matrix is computed iteratively stopping if the sum of the column sums of the current impact matrix increases less than a given threshold (colSumsThreshold) in comparison to the previous impact matrix. The resulting final impact matrix is saved using pre-defined file naming conventions. The resulting final impact matrix can be used to determine absolute impacts of source on target genes. Absolute impacts allow to identify those genes that have the greatest impact on other genes. Absolute impacts do not distinguish between potential activator and inhibitory contributions along the network paths. This function allows to separate the computation of the initial from the final impact matrix.
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation of impact matrix if improvement is less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}
#' @return Saves final impact matrix. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' computeNetworkFlow_CohortSpecificAbsoluteImpacts( dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath )
#' 
computeNetworkFlow_CohortSpecificAbsoluteImpacts <- function( dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE, patientSpecific = FALSE )
{
 
  if( output )
  {
      print( "Load basic network flow matrix:" )
  }
 
  ##Load network flow matrix
  subPath = "/NetworkPropagation/NetworkFlow/BasicNetworkFlowMatrices/"
  if( !patientSpecific )
  {
      basicNetworkFlowMatrixFile = paste( path, subPath, dataSetName, "_BasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".Rout", sep = "" )
  }
  else
  {
      basicNetworkFlowMatrixFile = paste( path, subPath, dataSetName, "_BasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".Rout", sep = "" )
  }
  
  #Loads objects: directCopyNumberContribution, networkFlowMatrix, ...
  load( file = basicNetworkFlowMatrixFile )
  
  ##Gives an error
  #print( colSums( as.matrix( networkFlowMatrix ) ) )  
  ##This helps, but would it be better to multiply a sparse matrix???
  networkFlowMatrix = as.matrix( networkFlowMatrix )
  
  ##Does not help, colSums gives still an error
  #networkFlowMatrix = Matrix( networkFlowMatrix )
  #networkFlowMatrix <- as( networkFlowMatrix, "sparseMatrix" )
  
  if( output )
  {
      print( basicNetworkFlowMatrixFile )
  }
  
  
  ##Compute network flow
  networkFlow = networkFlowMatrix
  pathwayContributions = networkFlowMatrix
  curColSums  = colSums( networkFlow )
  
  k = 1
  
  if( output )
  {
      print( paste( "Max column sum network flow matrix: ", max( curColSums ) ) )
      print( paste( "Min column norm network flow matrix: ", min( curColSums ) ) )
      print( paste( "CurColSums( ", k, " ): ", sum( curColSums ), sep = "" ) )
  }
  
  if( output )
  {
      print( "Compute network flow matrix" )
  }
  
  trackColSums = c( curColSums )
  trackColSumsDiff = c()
  repeat
  { 
    k = k + 1
    
    #Previous column sums to determine convergence of the network flow
    preColSums  = curColSums
    
    ##Compute contributions of pathways of length k     
    pathwayContributions = pathwayContributions %*% networkFlowMatrix 
    
    ##Add contributions of pathways of length k to network flow
    networkFlow = networkFlow + pathwayContributions
    
    #Current column sums
    curColSums  = colSums( networkFlow )
    trackColSums = c( trackColSums, curColSums )
    
    if( output )
    {
      print( paste( "CurColSums( ", k, " ): ", sum( curColSums ), sep = "" ) )
    }
    
    ##Convergence?
    colSumsDiff = sum( curColSums - preColSums )
    trackColSumsDiff = c( trackColSumsDiff, colSumsDiff )
    
    if( output )
    {
      print( paste( "ColSumsDiff( paths <= ", k, " ): ", colSumsDiff, sep = "" ) )
    }
    if( colSumsDiff < colSumsThreshold )
    {
      if( output )
      {
        print( paste( "Convergence reached after ", k, " steps.", sep = "" ) )
      }
      break
    }
  }
  
  if( output )
  {
      print( "Save network flow matrix:" )
  }
  
  ##Save output
  subPath = "/NetworkPropagation/NetworkFlow/FinalNetworkFlowMatrices/"
  if( !patientSpecific )
  {
      outputFile = paste( path, subPath, dataSetName, "_NetworkFlowMatrix_CohortSpecificAbsoluteImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, "_FlowCutoff_", colSumsThreshold, ".Rout", sep = "" )
  }
  else
  {
      outputFile = paste( path, subPath, dataSetName, "_NetworkFlowMatrix_PatientSpecificAbsoluteImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, "_FlowCutoff_", colSumsThreshold, ".Rout", sep = "" )
  }
  
  save( k, trackColSums, trackColSumsDiff, networkFlow, file = outputFile )
  
  if( output )
  {
      print( outputFile )
  }
  
}



#
#Load network flow matrix: internal function.
#
loadNetworkFlowMatrix_CohortSpecificAbsoluteImpacts <- function( dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = FALSE )
{
  subPath = "/NetworkPropagation/NetworkFlow/FinalNetworkFlowMatrices/"
  loadFile = paste( path, subPath, dataSetName, "_NetworkFlowMatrix_CohortSpecificAbsoluteImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, "_FlowCutoff_", colSumsThreshold, ".Rout", sep = "" )
  
  #Loads network flow objects: trackColSums, trackColSumsDiff, networkFlow
  load( loadFile )
  
  if( output )
  {
      print( loadFile )
  }
  
  return( networkFlow )
}



##
##
##
##
##
##
##
########################################################################################################################################################
########################################################################################################################################################
##
## Compute patient-specific absolute impacts
##
########################################################################################################################################################
########################################################################################################################################################
##
##
##
##
##
##
##

#' Network propagation module: Compute absolute patient-specific impact matrix
#'
#' This function computes an absolute patient-specific impact matrix for a patient (patient) part of a given data set (data) under consideration of a given network (networkName) using network propagation. Only significant regulatory links (pValCutoff) that are not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene are included in the computations. First, a basic initial patient-specific matrix is computed. Next, this matrix is used to iteratively compute the absolute patient-specific impact matrix stopping if the sum of the column sums of the current impact matrix increases less than a given threshold (colSumsThreshold) in comparison to the previous impact matrix. The basic inital impact matrix and the resulting final impact matrix are saved using pre-defined file naming conventions. The resulting final impact matrix can be used to determine absolute impacts of source on target genes. Absolute impacts allow to identify those genes that have the greatest impact on other genes. Absolute impacts do not distinguish between potential activator and inhibitory contributions along the network paths.
#' @param patient Patient considered for impact matrix computation
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation of impact matrix if improvement is less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlow_PatientSpecificAbsoluteImpacts}}, \code{\link{loadNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}, \code{\link{getAverageImpacts_PatientSpecificAbsoluteImpacts}}, \code{\link{getImpacts_PatientSpecificAbsoluteImpacts}}, \code{\link{getRawImpacts_PatientSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}}, ...
#' @return Saves basic initial and resulting final impact matrices. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#'
#' #Compute impact matrix for first patient of the given data set
#' computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts( patient = 1, data = data, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath )
#' 
computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts <- function( patient, data, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
    ##Compute basic network flow matrix
    computeBasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts( patient = patient, data = data, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, path = path, output = output )
    
    ##Compute final network flow matrix
    computeNetworkFlow_PatientSpecificAbsoluteImpacts( patient = patient, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
}



#' Network propagation module: Compute basic initial absolute patient-specific impact matrix
#'
#' This function computes a basic initial patient-specific impact matrix that is required for the iterative computation of an absolute patient-specific impact matrix. The basic matrix is computed for a given data set (data) under consideration of a given network (networkName) using network propagation. Only significant regulatory links (pValCutoff) that are not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene are included in the computations. The basic inital impact matrix is saved using pre-defined file naming conventions. This matrix can be used to compute an absolute patient-specific impact matrix using the function \code{\link{computeNetworkFlow_PatientSpecificAbsoluteImpacts}}. This function allows to separate the computation of the initial from the final impact matrix.
#' @param patient Patient considered for impact matrix computation
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlow_PatientSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}.
#' @return Saves basic initial impact matrix. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#'
#' #Compute basic initial impact matrix for first patient of the given data set
#' computeBasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts( patient = 1, data = data, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#' 
computeBasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts <- function( patient, data, dataSetName, networkName, pValCutoff, localGeneCutoff, path, output = TRUE )
{
  ##Load network and remove non-significant predictors
  network = loadNetworkWithFilteringForSignificantPredictors( networkName = networkName, pValCutoff = pValCutoff, path = path, output = output )
  
  ##Check if genes in the data set considered for prediction are in the same order of genes as in the data set used to learn the network
  if( !identical( data$genes, network$responseGenes ) )
  {
      print( paste( "Error: Genes in the given data set '", dataSetName, "' are not in the same order as the genes in the data set used to learn the network '", networkName, "'.", sep = "" ) )
      print( "" )
      print( " - Please ensure that the given data set has exactly the same order of genes." )
      print( " - Remove genes from the given data set if they were not contained in the data set used for network learning." )
      print( " - Add missing genes with artifical measurement values of zero to the given data set." )
      
      return()
  }  
  
  ##
  ##Basic variables
  ##    
  #Number of response variables: all genes in the data set
  nrResponses = length( data$genes )

  #Number of predictor variables (gene-specific copy number + expression levels)
  nrPredictors = nrow( data$Y ) + 1

  #Number of experiments
  T = ncol( data$U )  
  
  ##
  ##Basic matrices: Y and U contain data of considered response genes; M must be complete
  ##   
  #Response data
  #Gene expression levels of response genes: Y matrix -> number of experiments x number of response genes
  Y   = t( data.matrix( data$Y ) )
  
  #Predictor data
  #Copy numbers of response genes: U matrix -> number of experiments x number of response genes
  U   = t( data.matrix( data$U ) )      
  #Gene expression levels of predictor genes: M matrix -> number of experiments x number of predictor genes
  M   = t( data.matrix( data$Y ) )
  
  
 
  ##Network flow matrix
  #- rows: predictor genes
  #- colums: response genes
  #- rows and columns are exchanged compared to the network matrix to enable useful matrix multiplications
  #- square matrix: nrResponses x nrResponses to enable network flow analyses via matrix multiplication
  #
  
  networkFlowMatrix = matrix( 0, nrow = nrResponses, ncol = nrResponses )
  colnames( networkFlowMatrix ) <- data$genes
  rownames( networkFlowMatrix ) <- data$genes

  ##Direct copy number contributions
  directCopyNumberContribution = matrix( 0, nrow = nrResponses, ncol = 1 )
  rownames( directCopyNumberContribution ) <- data$genes
  colnames( directCopyNumberContribution ) <- "CopyNumberContribution"
  
  ##Explained variance
  correlationStatistics = matrix( 0, nrow = nrResponses, ncol = 1 )
  rownames( correlationStatistics ) <- data$genes
  colnames( correlationStatistics ) <- "AverageCorrelation"
  

  ##Explained variance
  explainedVariance = matrix( 0, nrow = nrResponses, ncol = 1 )
  rownames( explainedVariance ) <- data$genes
  colnames( explainedVariance ) <- "ExplainedVariance"

  if( output )
  {
      print( "Compute patient-specific basic network flow matrix" )
  }

  
  for( i in 1:nrResponses )
  {
      ##
      ##Response vector
      ##
      y = Y[ , i ]

      ##
      ##Remove predictors that are on the same chromosome like respone gene i in distance <= localGeneCutoff
      ##
      Mdummy = M
      
      if( localGeneCutoff > 0 )
      {
        currentResponseGene    = data$genes[ i ]
        currentResponseGeneChr = data$chr[ i ]

        #Determine downstream predictors
        removePredictorPos = c()
        for( v in max( 1, i - localGeneCutoff ):max( 1, i - 1 ) )
        {
          localNeighborGeneChr = data$chr[ v ]
   
          #If local neighbor is still on same chromosome then mark it for later removement
          if( localNeighborGeneChr == currentResponseGeneChr )
          {
            removePredictorPos = c( removePredictorPos, v )
          }      
        }

        #Determine upstream predictors
        for( v in min( nrResponses, i + 1 ):min( nrResponses, i + localGeneCutoff ) )
        {
          localNeighborGeneChr = data$chr[ v ]
  
          #If local neighbor is still on same chromosome then mark it for later removement
          if( localNeighborGeneChr == currentResponseGeneChr )
          {
            removePredictorPos = c( removePredictorPos, v )
          }      
        }

        #Remove predictors (note that the response gene itself is already removed during the network computation -> all network weights are zero for the response gene)
        Mdummy[ , removePredictorPos ] = 0      
      }
      
      ##
      ##Predictor matrix
      ##
      #Add copy number of response gene as a predictor
      specificM = cbind( U[ , i ], Mdummy )
                
      
      ##
      ##Compute contributions of individual predictors on the expression level of response gene i (contributions are patient-specific)
      ##      
      #Contributions of individual predictors (first entry is the copy number effect the remaining entries reflect effects from trans predictors)
      predictorContributions = computePredictorProportions_PatientSpecificAbsoluteImpacts( data = specificM[ patient, ], regressionParameters = network$G[ i,  ], predictorGenes = network$predictors, responseGenes = network$responseGenes )
    
      ##
      ##Compute the correlation between predicted and observed gene expression levels of gene i (correlations are cohort-specific)
      ##
      #Predict expression level ignoring offset (later correlation evaluation is invariant to the offset)
      Yhat = specificM %*%  network$G[ i,  ]
      
      corStatistics = 0
      if( var( Yhat, na.rm = TRUE ) == 0 )
      {
        if( output )
        {
          print( paste( data$genes[ i ], "not predictable" ) )
        }
      }
      else
      {      
        #Compute correlation and p-value
        dummy = cor.test( y, Yhat, alternative = "greater" )

        corStatistics = dummy[[ 4 ]] #Correlation
      }
          
      
      ##
      ##Compute the proportion of variance for gene i that is explained by the network
      ##            
      #Prediction not possible: e.g. no predictors available
      if( is.na( corStatistics ) )
      {
        corStatistics = 0
      }
      
      #Network cannot predict the expression level correctly
      if( corStatistics < 0 )
      {
        corStatistics = 0 
      }
      
      #Proportion of variance of gene i explained by the network
      squareCorrelation = corStatistics^2
      
      ##Multiply the predictor contributions with the explained variance
      predictorContributions = predictorContributions * squareCorrelation
      
      ##Transfer predictor contributions to the corresponding network flow column
      networkFlowMatrix[ , i ] = predictorContributions[ 2:( nrResponses + 1 ) ]

      ##
      ##Integrate direct copy number effect
      ##
      #- Entry [i,i] is always zero, since no gene can have a trans-contribution on its own expression
      #- We can use [i,i] to integrate the gene-specific copy number contribution 
      networkFlowMatrix[ i, i ] = predictorContributions[ 1 ]

      ##Transfer direct copy number contribution
      directCopyNumberContribution[ i, 1 ] = predictorContributions[ 1 ]
      
      ##Transfer average correlation
      correlationStatistics[ i, 1 ] = corStatistics
      
      ##Transfer explained variance
      explainedVariance[ i, 1 ] = squareCorrelation
  }

  ##
  ##Save basic network flow matrix
  ##
  if( output )
  {
      print( "Save output:" )
  }
  
  subPath = "/NetworkPropagation/NetworkFlow/BasicNetworkFlowMatrices/"
  outputFile =  paste( path, subPath, dataSetName, "_Patient_", patient, "_BasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".Rout", sep = "" )
  
  dummy <- Matrix( networkFlowMatrix, sparse = TRUE )
  networkFlowMatrix <- dummy
  
  save( networkFlowMatrix, directCopyNumberContribution, explainedVariance, correlationStatistics, file = outputFile )
  
  if( output )
  {
      print( outputFile )
  }
}



#
#Compute the proportion that individual predictors have on the prediction of the expression level of a gene in a specific patient: internal function.
#
# data: patient-specific vector (copy number of considered gene 'i' above and expression levels of all other genes)
#
computePredictorProportions_PatientSpecificAbsoluteImpacts <- function( data, regressionParameters, predictorGenes, responseGenes )
{
    activePredictors = which( regressionParameters != 0 )    
    N = length( activePredictors )
    
    predictorProportions = rep( 0, N )
    
    #Global return values
    copyNumberContribution = 0 
    globalTransPredictorContributionVector = rep( 0, times = length( responseGenes ) )
    
    ##
    ##Compute absolute proportion of individual predictors
    ##
    if( N > 0 )
    {
      ##Compute proportion for the specific 'patient'
      for( i in 1:N )
      {
        #Average absolute contribution of the predictor is obtained by
        #- multiplying the measurements of the predictor with its corresponding regression parameter
        #- followed by averaging over the different tumor samples
        predictorProportions[ i ] = abs( data[ activePredictors[ i ] ] * regressionParameters[ activePredictors[ i ] ] )	  
      }
      
      #Standardization
      predictorProportions = predictorProportions / sum( predictorProportions, na.rm = TRUE )
      
      ##
      ##Create global predictor contribution vector
      ##
      #Determine copy number contribution if active and remove copy number predictor from the list of active predictors, because we are later interested in the influence of trans effects         
      if( activePredictors[ 1 ] == 1 )
      {
        copyNumberContribution = predictorProportions[ 1 ]  
        activePredictors = setdiff( activePredictors, 1 )
        if( N > 1 )
        {
          #remove copy number proportion
          predictorProportions = predictorProportions[ 2:N ]
        }
      }
      
      ##This shift is necessary, since the regression parameters start with the copy number predictor
      activePredictors = activePredictors - 1      
      
      N = length( activePredictors ) #if there is only a copy number predictor, N will be zero
      if( N > 0 )
      {
        activePredictorGenes = responseGenes[ activePredictors ]  
        mapPos = match( activePredictorGenes, responseGenes )
  
        ##Map predictor genes to their corresponding position in the trans-contribution vector (order: responseGenes)
        uniqueMapPos = mapPos[ which( !is.na( mapPos ) ) ]
  
        globalTransPredictorContributionVector[ uniqueMapPos ] = predictorProportions[ which( !is.na( mapPos ) ) ]	    
      }      
      
    }
    
    res  = checkContributions( res = c( copyNumberContribution, globalTransPredictorContributionVector ) )
   
    #First entry is the copy number contribution the remaining entries reflect the trans predictors
    return( res )
} 



#' Network propagation module: Compute absolute patient-specific impact matrix
#'
#' This function computes an absolute patient-specific impact matrix for a patient (patient) in a given data set (dataSetName) under consideration of a given network (networkName) using network propagation. The basis of this computation is the pre-computed basic initial patient-specific impact matrix that only considered significant regulatory links (pValCutoff) that were not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene (\code{\link{computeBasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}). The absolute patient-specific impact matrix is computed iteratively stopping if the sum of the column sums of the current impact matrix increases less than a given threshold (colSumsThreshold) in comparison to the previous impact matrix. The resulting final impact matrix is saved using pre-defined file naming conventions. The resulting final impact matrix can be used to determine absolute impacts of source on target genes. Absolute impacts allow to identify those genes that have the greatest impact on other genes. Absolute impacts do not distinguish between potential activator and inhibitory contributions along the network paths. This function allows to separate the computation of the initial from the final impact matrix.
#' @param patient Patient considered for impact matrix computation
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation of impact matrix if improvement is less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}
#' @return Saves final impact matrix. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' computeNetworkFlow_PatientSpecificAbsoluteImpacts( patient = 1, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath )
#'
computeNetworkFlow_PatientSpecificAbsoluteImpacts <- function( patient, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
  ##Network flow computation is identical with cohort-specific computations, but now the patient-specific basic network flow matrix is utilized
  computeNetworkFlow_CohortSpecificAbsoluteImpacts( dataSetName = paste( dataSetName, "_Patient_", patient, sep = "" ), networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output, patientSpecific = TRUE )
}



#
#Load network flow matrix: internal function.
#
loadNetworkFlowMatrix_PatientSpecificAbsoluteImpacts <- function( patient, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = FALSE )
{
  subPath = "/NetworkPropagation/NetworkFlow/FinalNetworkFlowMatrices/"
  loadFile = paste( path, subPath, dataSetName, "_Patient_", patient, "_NetworkFlowMatrix_PatientSpecificAbsoluteImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, "_FlowCutoff_", colSumsThreshold, ".Rout", sep = "" )
  
  #Loads network flow objects: trackColSums, trackColSumsDiff, networkFlow
  load( loadFile )

  if( output )
  {
      print( loadFile )
  }

  
  return( networkFlow )
}



##
##
##
##
##
##
##
########################################################################################################################################################
########################################################################################################################################################
##
## Compute patient-specific relative impacts
##
########################################################################################################################################################
########################################################################################################################################################
##
##
##
##
##
##
##

#' Network propagation module: Compute relative patient-specific impact matrix
#'
#' This function computes a relative patient-specific impact matrix for a patient (patient) part of a given data set (data) under consideration of a given network (networkName) using network propagation. Only significant regulatory links (pValCutoff) that are not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene are included in the computations. First, a basic initial patient-specific matrix is computed. Next, this matrix is used to iteratively compute the relative patient-specific impact matrix stopping if the sum of the column sums of the current impact matrix and the previous impact matrix differ less than a given threshold (colSumsThreshold). The basic inital impact matrix and the resulting final impact matrix are saved using pre-defined file naming conventions. The resulting final impact matrix can be used to determine relative impacts of source on target genes. Relative impacts account for potential activator and inhibitory contributions along the network paths.
#' @param patient Patient considered for impact matrix computation
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation of impact matrix if global change of impacts is less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlow_PatientSpecificAbsoluteImpacts}}, \code{\link{loadNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}, \code{\link{getAverageImpacts_PatientSpecificAbsoluteImpacts}}, \code{\link{getImpacts_PatientSpecificAbsoluteImpacts}}, \code{\link{getRawImpacts_PatientSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}}, ...
#' @return Saves basic initial and resulting final impact matrices. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#'
#' #Compute impact matrix for first patient of the given data set
#' computeNetworkFlowMatrix_PatientSpecificRelativeImpacts( patient = 1, data = data, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath )
#' 
computeNetworkFlowMatrix_PatientSpecificRelativeImpacts <- function( patient, data, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
    ##Compute basic network flow matrix
    computeBasicNetworkFlowMatrix_PatientSpecificRelativeImpacts( patient = patient, data = data, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, path = path, output = output )
    
    ##Compute final network flow matrix
    computeNetworkFlow_PatientSpecificRelativeImpacts( patient = patient, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
}



#' Network propagation module: Compute basic initial relative patient-specific impact matrix
#'
#' This function computes a basic initial patient-specific impact matrix that is required for the iterative computation of a relative patient-specific impact matrix. The basic matrix is computed for a given data set (data) under consideration of a given network (networkName) using network propagation. Only significant regulatory links (pValCutoff) that are not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene are included in the computations. The basic inital impact matrix is saved using pre-defined file naming conventions. This matrix can be used to compute a relative patient-specific impact matrix using the function \code{\link{computeNetworkFlow_PatientSpecificRelativeImpacts}}. This function allows to separate the computation of the initial from the final impact matrix.
#' @param patient Patient considered for impact matrix computation
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlow_PatientSpecificRelativeImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}}.
#' @return Saves basic initial impact matrix. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#'
#' #Compute basic initial impact matrix for first patient of the given data set
#' computeBasicNetworkFlowMatrix_PatientSpecificRelativeImpacts( patient = 1, data = data, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#' 
computeBasicNetworkFlowMatrix_PatientSpecificRelativeImpacts <- function( patient, data, dataSetName, networkName, pValCutoff, localGeneCutoff, path, output = TRUE )
{
  ##Load network and remove non-significant predictors
  network = loadNetworkWithFilteringForSignificantPredictors( networkName = networkName, pValCutoff = pValCutoff, path = path, output = output )

  ##Check if genes in the data set considered for prediction are in the same order of genes as in the data set used to learn the network
  if( !identical( data$genes, network$responseGenes ) )
  {
      print( paste( "Error: Genes in the given data set '", dataSetName, "' are not in the same order as the genes in the data set used to learn the network '", networkName, "'.", sep = "" ) )
      print( "" )
      print( " - Please ensure that the given data set has exactly the same order of genes." )
      print( " - Remove genes from the given data set if they were not contained in the data set used for network learning." )
      print( " - Add missing genes with artifical measurement values of zero to the given data set." )
      
      return()
  }  
  
  ##
  ##Basic variables
  ##    
  #Number of response variables: all genes in the data set
  nrResponses = length( data$genes )

  #Number of predictor variables (gene-specific copy number + expression levels)
  nrPredictors = nrow( data$Y ) + 1

  #Number of experiments
  T = ncol( data$U )  
  
  ##
  ##Basic matrices: Y and U contain data of considered response genes; M must be complete
  ##   
  #Response data
  #Gene expression levels of response genes: Y matrix -> number of experiments x number of response genes
  Y   = t( data.matrix( data$Y ) )
  
  #Predictor data
  #Copy numbers of response genes: U matrix -> number of experiments x number of response genes
  U   = t( data.matrix( data$U ) )      
  #Gene expression levels of predictor genes: M matrix -> number of experiments x number of predictor genes
  M   = t( data.matrix( data$Y ) )
  
  
 
  ##Network flow matrix
  #- rows: predictor genes
  #- colums: response genes
  #- rows and columns are exchanged compared to the network matrix to enable useful matrix multiplications
  #- square matrix: nrResponses x nrResponses to enable network flow analyses via matrix multiplication
  #
  
  networkFlowMatrix = matrix( 0, nrow = nrResponses, ncol = nrResponses )
  colnames( networkFlowMatrix ) <- data$genes
  rownames( networkFlowMatrix ) <- data$genes

  ##Direct copy number contributions
  directCopyNumberContribution = matrix( 0, nrow = nrResponses, ncol = 1 )
  rownames( directCopyNumberContribution ) <- data$genes
  colnames( directCopyNumberContribution ) <- "CopyNumberContribution"
  
  ##Explained variance
  correlationStatistics = matrix( 0, nrow = nrResponses, ncol = 1 )
  rownames( correlationStatistics ) <- data$genes
  colnames( correlationStatistics ) <- "AverageCorrelation"
  

  ##Explained variance
  explainedVariance = matrix( 0, nrow = nrResponses, ncol = 1 )
  rownames( explainedVariance ) <- data$genes
  colnames( explainedVariance ) <- "ExplainedVariance"

  if( output )
  {
      print( "Compute patient-specific basic network flow matrix" )
  }

  
  for( i in 1:nrResponses )
  {
      ##
      ##Response vector
      ##
      y = Y[ , i ]

      ##
      ##Remove predictors that are on the same chromosome like respone gene i in distance <= localGeneCutoff
      ##
      Mdummy = M
      
      if( localGeneCutoff > 0 )
      {
        currentResponseGene    = data$genes[ i ]
        currentResponseGeneChr = data$chr[ i ]

        #Determine downstream predictors
        removePredictorPos = c()
        for( v in max( 1, i - localGeneCutoff ):max( 1, i - 1 ) )
        {
          localNeighborGeneChr = data$chr[ v ]
  
          #If local neighbor is still on same chromosome then mark it for later removement
          if( localNeighborGeneChr == currentResponseGeneChr )
          {
            removePredictorPos = c( removePredictorPos, v )
          }      
        }

        #Determine upstream predictors
        for( v in min( nrResponses, i + 1 ):min( nrResponses, i + localGeneCutoff ) )
        {
          localNeighborGeneChr = data$chr[ v ]
    
          #If local neighbor is still on same chromosome then mark it for later removement
          if( localNeighborGeneChr == currentResponseGeneChr )
          {
            removePredictorPos = c( removePredictorPos, v )
          }      
        }

        #Remove predictors (note that the response gene itself is already removed during the network computation -> all network weights are zero for the response gene)
        Mdummy[ , removePredictorPos ] = 0      
      }
      
      
      ##
      ##Predictor matrix
      ##
      #Add copy number of response gene as a predictor
      specificM = cbind( U[ , i ], Mdummy )
                
      
      ##
      ##Compute contributions of individual predictors on the expression level of response gene i (contributions are patient-specific and distinguish between repression and activation)
      ##      
      #Contributions of individual predictors (first entry is the copy number effect the remaining entries reflect effects from trans predictors)
      predictorContributions = computePredictorProportions_PatientSpecificRelativeImpacts( data = specificM[ patient, ], regressionParameters = network$G[ i,  ], predictorGenes = network$predictors, responseGenes = network$responseGenes )
    
      ##
      ##Compute the correlation between predicted and observed gene expression levels of gene i (correlations are cohort-specific)
      ##
      #Predict expression level ignoring offset (later correlation evaluation is invariant to the offset)
      Yhat = specificM %*%  network$G[ i,  ]
      
      corStatistics = 0
      if( var( Yhat, na.rm = TRUE ) == 0 )
      {
        if( output )
        {
          print( paste( data$genes[ i ], "not predictable" ) )
        }
      }
      else
      {      
        #Compute correlation and p-value
        dummy = cor.test( y, Yhat, alternative = "greater" )

        corStatistics = dummy[[ 4 ]] #Correlation
      }
          
      
      ##
      ##Compute the proportion of variance for gene i that is explained by the network
      ##            
      #Prediction not possible: e.g. no predictors available
      if( is.na( corStatistics ) )
      {
	corStatistics = 0
      }
      
      #Network cannot predict the expression level correctly
      if( corStatistics < 0 )
      {
	  corStatistics = 0 
      }
      
      #Proportion of variance of gene i explained by the network
      squareCorrelation = corStatistics^2
      
      ##Multiply the predictor contributions with the explained variance
      predictorContributions = predictorContributions * squareCorrelation
      
      ##Transfer predictor contributions to the corresponding network flow column
      networkFlowMatrix[ , i ] = predictorContributions[ 2:( nrResponses + 1 ) ]

      ##
      ##Integrate direct copy number effect
      ##
      #- Entry [i,i] is always zero, since no gene can have a trans-contribution on its own expression
      #- We can use [i,i] to integrate the gene-specific copy number contribution 
      networkFlowMatrix[ i, i ] = predictorContributions[ 1 ]

      ##Transfer direct copy number contribution
      directCopyNumberContribution[ i, 1 ] = predictorContributions[ 1 ]
      
      ##Transfer average correlation
      correlationStatistics[ i, 1 ] = corStatistics
      
      ##Transfer explained variance
      explainedVariance[ i, 1 ] = squareCorrelation
  }

  ##
  ##Save basic network flow matrix
  ##
  if( output )
  {
      print( "Save output:" )
  }
  
  subPath = "/NetworkPropagation/NetworkFlow/BasicNetworkFlowMatrices/"
  outputFile =  paste( path, subPath, dataSetName, "_Patient_", patient, "_BasicNetworkFlowMatrix_PatientSpecificRelativeImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".Rout", sep = "" )
  
  dummy <- Matrix( networkFlowMatrix, sparse = TRUE )
  networkFlowMatrix <- dummy
  
  save( networkFlowMatrix, directCopyNumberContribution, explainedVariance, correlationStatistics, file = outputFile )
  
  if( output )
  {
      print( outputFile )
  }

}


#
#Compute the proportion that individual predictors have on the prediction of the expression level of a gene in a specific patient: internal function.
#
computePredictorProportions_PatientSpecificRelativeImpacts <- function( data, regressionParameters, predictorGenes, responseGenes )
{
    activePredictors = which( regressionParameters != 0 )    
    N = length( activePredictors )
    
    predictorProportions = rep( 0, N )
    predictorProportionsWithSign = rep( 0, N )  

    #Global return values
    copyNumberContribution = 0 
    globalTransPredictorContributionVector = rep( 0, times = length( responseGenes ) )
    
    ##
    ##Compute absolute proportion of individual predictors
    ##
    if( N > 0 )
    {
      ##Compute proportion for the specific 'patient'
      for( i in 1:N )
      {
        #Average absolute contribution of the predictor is obtained by
        #- multiplying the measurements of the predictor with its corresponding regression parameter
        #- followed by averaging over the different tumor samples
        predictorProportions[ i ] = abs( data[ activePredictors[ i ] ] * regressionParameters[ activePredictors[ i ] ] )

        #Average contribution of the predictor with directionality (average is later determined via standardization)
        predictorProportionsWithSign[ i ] = data[ activePredictors[ i ] ] * regressionParameters[ activePredictors[ i ] ]  
      }
      
      #Standardization: Keep directionality of predictor effects due to the usage of predictorProportionsWithSign
      predictorProportions = predictorProportionsWithSign / sum( predictorProportions, na.rm = TRUE )
      
      ##
      ##Create global predictor contribution vector
      ##
      #Determine copy number contribution if active and remove copy number predictor from the list of active predictors, because we are later interested in the influence of trans effects         
      if( activePredictors[ 1 ] == 1 )
      {
        copyNumberContribution = predictorProportions[ 1 ]  
        activePredictors = setdiff( activePredictors, 1 )
        if( N > 1 )
        {
          #remove copy number proportion
          predictorProportions = predictorProportions[ 2:N ]
        }
      }
      
      ##This shift is necessary, since the regression parameters start with the copy number predictor
      activePredictors = activePredictors - 1      
      
      N = length( activePredictors ) #if there is only a copy number predictor, N will be zero
      if( N > 0 )
      {
        activePredictorGenes = responseGenes[ activePredictors ]  
        mapPos = match( activePredictorGenes, responseGenes )
  
        ##Map predictor genes to their corresponding position in the trans-contribution vector (order: responseGenes)
        uniqueMapPos = mapPos[ which( !is.na( mapPos ) ) ]
 
        globalTransPredictorContributionVector[ uniqueMapPos ] = predictorProportions[ which( !is.na( mapPos ) ) ]  
      }      
      
    }
    
    res  = checkContributions( res = c( copyNumberContribution, globalTransPredictorContributionVector ) )
   
    #First entry is the copy number contribution the remaining entries reflect the trans predictors
    return( res )
} 



#' Network propagation module: Compute relative patient-specific impact matrix
#'
#' This function computes a relative patient-specific impact matrix for a patient (patient) in a given data set (dataSetName) under consideration of a given network (networkName) using network propagation. The basis of this computation is the pre-computed basic initial patient-specific impact matrix that only considered significant regulatory links (pValCutoff) that were not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene (\code{\link{computeBasicNetworkFlowMatrix_PatientSpecificRelativeImpacts}}). The relative patient-specific impact matrix is computed iteratively stopping if the sum of the column sums of the current impact matrix and the previous impact matrix differ less than a given threshold (colSumsThreshold). The resulting final impact matrix is saved using pre-defined file naming conventions. The resulting final impact matrix can be used to determine relative impacts of source on target genes. Relative impacts account for potential activator and inhibitory contributions along the network paths. This function allows to separate the computation of the initial from the final impact matrix.
#' @param patient Patient considered for impact matrix computation
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation of impact matrix if global change of impacts is less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}
#' @return Saves final impact matrix. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' computeNetworkFlow_PatientSpecificRelativeImpacts( patient = 1, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath )
#'
computeNetworkFlow_PatientSpecificRelativeImpacts <- function( patient, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
 
  if( output )
  {
    print( "Load basic network flow matrix:" )
  }
 
  subPath = "/NetworkPropagation/NetworkFlow/BasicNetworkFlowMatrices/"
  basicNetworkFlowMatrixFile = paste( path, subPath, dataSetName, "_Patient_", patient, "_BasicNetworkFlowMatrix_PatientSpecificRelativeImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".Rout", sep = "" )
  
  #Loads objects: directCopyNumberContribution, networkFlowMatrix, ...
  load( file = basicNetworkFlowMatrixFile )
  
  networkFlowMatrix = as.matrix( networkFlowMatrix )
  
  if( output )
  {
    print( basicNetworkFlowMatrixFile )
  }
  
  
  ##Compute network flow
  networkFlow = networkFlowMatrix
  pathwayContributions = networkFlowMatrix
  curColSums  = colSums( networkFlow )
  
  k = 1
  
  if( output )
  {
      print( paste( "Max column sum network flow matrix: ", max( curColSums ) ) )
      print( paste( "Min column norm network flow matrix: ", min( curColSums ) ) )
      print( paste( "CurColSums( ", k, " ): ", sum( curColSums ), sep = "" ) )
  }
  
  if( output )
  {
      print( "Compute network flow matrix" )
  }
  
  trackColSums = c( curColSums )
  trackColSumsDiff = c()
  repeat
  { 
    k = k + 1
    
    #Previous column sums to determine convergence of the network flow
    preColSums  = curColSums
    
    ##Compute contributions of pathways of length k     
    pathwayContributions = pathwayContributions %*% networkFlowMatrix 
    
    ##Add contributions of pathways of length k to network flow
    networkFlow = networkFlow + pathwayContributions
    
    #Current column sums
    curColSums  = colSums( networkFlow )
    trackColSums = c( trackColSums, curColSums )
    if( output )
    {
      print( paste( "CurColSums( ", k, " ): ", sum( curColSums ), sep = "" ) )
    }
    
    ##Convergence?
    colSumsDiff = abs( sum( curColSums - preColSums ) )
    trackColSumsDiff = c( trackColSumsDiff, colSumsDiff )
    if( output )
    {
      print( paste( "AbsoluteColSumsDiff( paths <= ", k, " ): ", colSumsDiff, sep = "" ) )
    }
    if( colSumsDiff < colSumsThreshold )
    {
      if( output )
      {
        print( paste( "Convergence reached after ", k, " steps.", sep = "" ) )
      }
      break
    }
  }
  
  if( output )
  {
      print( "Save network flow matrix:" )
  }
  
  ##Save output
  subPath = "/NetworkPropagation/NetworkFlow/FinalNetworkFlowMatrices/"
  outputFile = paste( path, subPath, dataSetName, "_Patient_", patient, "_NetworkFlowMatrix_PatientSpecificRelativeImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, "_FlowCutoff_", colSumsThreshold, ".Rout", sep = "" )
  
  save( k, trackColSums, trackColSumsDiff, networkFlow, file = outputFile )
  
  if( output )
  {
    print( outputFile )
  }
}



#
#Load network flow matrix: internal function.
#
loadNetworkFlowMatrix_PatientSpecificRelativeImpacts <- function( patient, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = FALSE )
{
  subPath = "/NetworkPropagation/NetworkFlow/FinalNetworkFlowMatrices/"
  loadFile = paste( path, subPath, dataSetName, "_Patient_", patient, "_NetworkFlowMatrix_PatientSpecificRelativeImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, "_FlowCutoff_", colSumsThreshold, ".Rout", sep = "" )
  
  #Loads network flow objects: trackColSums, trackColSumsDiff, networkFlow
  load( loadFile )

  if( output )
  {
      print( loadFile )
  }
  
  return( networkFlow )
}
