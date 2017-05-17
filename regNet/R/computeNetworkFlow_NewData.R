
##
##
##
##
########################################################################################################################################################
########################################################################################################################################################
##
## Compute cohort-specific absolute impacts for new data utilizing correlation statistics from a given cohort
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

#' Network propagation module: Compute basic initial absolute cohort-specific impact matrix
#'
#' This function computes a basic initial absolute cohort-specific impact matrix utilizing a given correlation statistics. This function provides the basis for the function \code{\link{computeNetworkFlowBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}}. This function allows to separate the computation of the initial from the final impact matrix. See \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}} and \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}} for more details.
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics are taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}
#' @return Saves basic initial impact matrix. See output for details.
#' @export
#' @examples
#'
#' #Load test and train data set and make test data set compatible to the train data set
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' gbmTFData = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#'
#' #Determine correlation statistics
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' predictGeneExpression( data = trainData, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
#' #Compute impact matrix using pre-computed correlation statistics
#' computeBasicNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( data = testData, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
computeBasicNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts <- function( data, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, path, output = TRUE )
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
  
  ##Correlation statistics (Use correlations between network-based predicted and originally measured expression levels from 'corStatDataSetName' that was already analyzed by the same network.)
  subPath = "/NetworkPredictions/"
  corStatDataFile = paste( path, subPath, corStatDataSetName, "_PredictionOfGeneExpressionBasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff ,".txt", sep = "" )
  corStat = read.delim( file = corStatDataFile, header = TRUE )#Gene, Correlation, P-Value
  correlationStatistics = matrix( as.numeric( corStat[ , 2 ] ), nrow = nrResponses, ncol = 1 )
  rownames( correlationStatistics ) <- as.character( corStat[ , 1 ] )
  colnames( correlationStatistics ) <- "AverageCorrelation"
  

  ##Explained variance
  explainedVariance = matrix( as.numeric( corStat[ , 2 ] )^2, nrow = nrResponses, ncol = 1 )
  rownames( explainedVariance ) <- as.character( corStat[ , 1 ] )
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
      ##Utilize the correlation between predicted and observed gene expression levels of gene i from 'corStatDataSetName'
      ##      
      corStatistics = correlationStatistics[ i, 1 ]
      if( is.na( corStatistics ) || corStatistics <= 0 )
      {
        if( output )
        {
          print( paste( data$genes[ i ], "not predictable" ) )
        }
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
  }
  

  ##
  ##Save basic network flow matrix
  ##
  if( output )
  {
      print( "Save output:" )
  }

  subPath = "/NetworkPropagation/NetworkFlow/BasicNetworkFlowMatrices/"
  outputFile =  paste( path, subPath, dataSetName, "_UsingCorStat_", corStatDataSetName, "_BasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".Rout", sep = "" )
  
  dummy <- Matrix( networkFlowMatrix, sparse = TRUE )
  networkFlowMatrix <- dummy
  
  save( networkFlowMatrix, directCopyNumberContribution, explainedVariance, correlationStatistics, file = outputFile )
  
  if( output )
  {
      print( outputFile )
  }

}



#' Network propagation module: Compute absolute cohort-specific impact matrix
#'
#' This function computes an absolute cohort-specific impact matrix utilizing a given correlation statistics. This function requires a pre-computed basic initial absolute cohort-specific impact matrix computed by the function \code{\link{computeBasicNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}}. This function allows to separate the computation of the initial from the final impact matrix. See \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}} and \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}} for more details.
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics are taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation of impact matrix if improvement is less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}
#' @return Saves final impact matrix. See output for details.
#' @export
#' @examples
#'
#' #Compute impact matrix using pre-computed correlation statistics
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' computeNetworkFlowBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
computeNetworkFlowBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts <- function( dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
  computeNetworkFlow_CohortSpecificAbsoluteImpacts( dataSetName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" ), networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output, patientSpecific = FALSE )
}



#' Network propagation module: Compute absolute cohort-specific impact matrix
#'
#' This function computes an absolute cohort-specific impact matrix utilizing a given correlation statistics. This enables to analyze a new data set (data, dataSetName) under consideration of the predictive power of another data set (corStatDataSetName). This supports the analysis of data sets with few patients where the sample size is too small to obtain a robust estimate of the predictive power of individual genes (correlation statistics). See \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}} for more details. The correlation statistics can be pre-computed using the function \code{\link{predictGeneExpression}}.
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics are taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation of impact matrix if improvement is less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}
#' @return Saves basic initial and resulting final impact matrices. See output for details.
#' @export
#' @examples
#'
#' #Load test and train data set and make test data set compatible to the train data set
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' gbmTFData = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#'
#' #Determine correlation statistics
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' predictGeneExpression( data = trainData, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
#' #Compute impact matrix using pre-computed correlation statistics
#' computeNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( data = testData, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath )
#'
computeNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts <- function( data, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
    ##Compute basic network flow matrix
    computeBasicNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( data = data, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, path = path, output = output )
    
    ##Compute final network flow matrix
    computeNetworkFlowBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
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
## Compute patient-specific absolute impacts for new data utilizing correlation statistics from a given cohort
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

#' Network propagation module: Compute basic initial absolute patient-specific impact matrix
#'
#' This function computes a basic initial absolute patient-specific impact matrix utilizing a given correlation statistics. This function provides the basis for the function \code{\link{computeNetworkFlowBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}}. This function allows to separate the computation of the initial from the final impact matrix. See \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}} and \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}} for more details.
#' @param patient Patient considered for impact matrix computation
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics are taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}
#' @return Saves basic initial impact matrix. See output for details.
#' @export
#' @examples
#'
#' #Load test and train data set and make test data set compatible to the train data set
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' gbmTFData = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#'
#' #Determine correlation statistics
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' predictGeneExpression( data = trainData, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
#' #Compute impact matrix for first patient using pre-computed correlation statistics
#' computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = 1, data = testData, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts <- function( patient, data, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, path, output = TRUE )
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
  
  ##Correlation statistics (Use correlations between network-based predicted and originally measured expression levels from 'corStatDataSetName' that was already analyzed by the same network.)
  subPath = "/NetworkPredictions/"
  corStatDataFile = paste( path, subPath, corStatDataSetName, "_PredictionOfGeneExpressionBasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff ,".txt", sep = "" )
  corStat = read.delim( file = corStatDataFile, header = TRUE )#Gene, Correlation, P-Value
  correlationStatistics = matrix( as.numeric( corStat[ , 2 ] ), nrow = nrResponses, ncol = 1 )
  rownames( correlationStatistics ) <- as.character( corStat[ , 1 ] )
  colnames( correlationStatistics ) <- "AverageCorrelation"
  

  ##Explained variance
  explainedVariance = matrix( as.numeric( corStat[ , 2 ] )^2, nrow = nrResponses, ncol = 1 )
  rownames( explainedVariance ) <- as.character( corStat[ , 1 ] )
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
      ##Utilize the correlation between predicted and observed gene expression levels of gene i from 'corStatDataSetName'
      ##      
      corStatistics = correlationStatistics[ i, 1 ]
      if( is.na( corStatistics ) || corStatistics <= 0 )
      {
        if( output )
        {
          print( paste( data$genes[ i ], "not predictable" ) )
        }
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
  }

  ##
  ##Save basic network flow matrix
  ##
  if( output )
  {
      print( "Save output:" )
  }
  
  subPath = "/NetworkPropagation/NetworkFlow/BasicNetworkFlowMatrices/"
  outputFile =  paste( path, subPath, dataSetName, "_UsingCorStat_", corStatDataSetName, "_Patient_", patient, "_BasicNetworkFlowMatrix_PatientSpecificAbsoluteImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".Rout", sep = "" )
  
  dummy <- Matrix( networkFlowMatrix, sparse = TRUE )
  networkFlowMatrix <- dummy
  
  save( networkFlowMatrix, directCopyNumberContribution, explainedVariance, correlationStatistics, file = outputFile )
  
  if( output )
  {
      print( outputFile )
  }

}



#' Network propagation module: Compute absolute patient-specific impact matrix
#'
#' This function computes an absolute patient-specific impact matrix utilizing a given correlation statistics. This function requires a pre-computed basic initial absolute patient-specific impact matrix computed by the function \code{\link{computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}}. This function allows to separate the computation of the initial from the final impact matrix. See \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}} and \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}} for more details.
#' @param patient Patient considered for impact matrix computation
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics are taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation of impact matrix if improvement is less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}
#' @return Saves final impact matrix. See output for details.
#' @export
#' @examples
#'
#' #Compute impact matrix for first patient using pre-computed correlation statistics
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' computeNetworkFlowBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = 1, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
computeNetworkFlowBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts <- function( patient, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
  ##Network flow computation is identical with cohort-specific computations, but now the patient-specific basic network flow matrix is utilized
  computeNetworkFlow_CohortSpecificAbsoluteImpacts( dataSetName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, "_Patient_", patient, sep = "" ), networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output, patientSpecific = TRUE )
}



#' Network propagation module: Compute absolute patient-specific impact matrix
#'
#' This function computes an absolute patient-specific impact matrix utilizing a given correlation statistics. This enables to analyze a patient (patient) in a new data set (data, dataSetName) under consideration of the predictive power of another data set (corStatDataSetName). This supports the analysis of data sets with few patients where the sample size is too small to obtain a robust estimate of the predictive power of individual genes (correlation statistics). See \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}} for more details. The correlation statistics can be pre-computed using the function \code{\link{predictGeneExpression}}.
#' @param patient Patient considered for impact matrix computation
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics are taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation of impact matrix if improvement is less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}
#' @return Saves basic initial and resulting final impact matrices. See output for details.
#' @export
#' @examples
#'
#' #Load test and train data set and make test data set compatible to the train data set
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' gbmTFData = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#'
#' #Determine correlation statistics
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' predictGeneExpression( data = trainData, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
#' #Compute impact matrix for first patient using pre-computed correlation statistics
#' computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = 1, data = testData, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath )
#'
computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts <- function( patient, data, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
    ##Compute basic network flow matrix
    computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = patient, data = data, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, path = path, output = output )
    
    ##Compute final network flow matrix
    computeNetworkFlowBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = patient, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
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
## Compute patient-specific relative impacts for new data utilizing correlation statistics from a given cohort
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

#' Network propagation module: Compute basic initial relative patient-specific impact matrix
#'
#' This function computes a basic initial relative patient-specific impact matrix utilizing a given correlation statistics. This function provides the basis for the function \code{\link{computeNetworkFlowBasedOnGivenCorStat_PatientSpecificRelativeImpacts}}. This function allows to separate the computation of the initial from the final impact matrix. See \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts}} and \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}} for more details.
#' @param patient Patient considered for impact matrix computation
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics are taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}}
#' @return Saves basic initial impact matrix. See output for details.
#' @export
#' @examples
#'
#' #Load test and train data set and make test data set compatible to the train data set
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' gbmTFData = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#'
#' #Determine correlation statistics
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' predictGeneExpression( data = trainData, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
#' #Compute impact matrix for first patient using pre-computed correlation statistics
#' computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = 1, data = testData, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts <- function( patient, data, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, path, output = TRUE )
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
  
  ##Correlation statistics (Use correlations between network-based predicted and originally measured expression levels from 'corStatDataSetName' that was already analyzed by the same network.)
  subPath = "/NetworkPredictions/"
  corStatDataFile = paste( path, subPath, corStatDataSetName, "_PredictionOfGeneExpressionBasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff ,".txt", sep = "" )
  corStat = read.delim( file = corStatDataFile, header = TRUE )#Gene, Correlation, P-Value
  correlationStatistics = matrix( as.numeric( corStat[ , 2 ] ), nrow = nrResponses, ncol = 1 )
  rownames( correlationStatistics ) <- as.character( corStat[ , 1 ] )
  colnames( correlationStatistics ) <- "AverageCorrelation"

  ##Explained variance
  explainedVariance = matrix( as.numeric( corStat[ , 2 ] )^2, nrow = nrResponses, ncol = 1 )
  rownames( explainedVariance ) <- as.character( corStat[ , 1 ] )
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
      ##Utilize the correlation between predicted and observed gene expression levels of gene i from 'corStatDataSetName'
      ##      
      corStatistics = correlationStatistics[ i, 1 ]
      if( is.na( corStatistics ) || corStatistics <= 0 )
      {
        if( output )
        {
          print( paste( data$genes[ i ], "not predictable" ) )
        }
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
  }

  ##
  ##Save basic network flow matrix
  ##
  if( output )
  {
      print( "Save output:" )
  }
  
  subPath = "/NetworkPropagation/NetworkFlow/BasicNetworkFlowMatrices/"
  outputFile =  paste( path, subPath, dataSetName, "_UsingCorStat_", corStatDataSetName, "_Patient_", patient, "_BasicNetworkFlowMatrix_PatientSpecificRelativeImpacts_BasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".Rout", sep = "" )
  
  dummy <- Matrix( networkFlowMatrix, sparse = TRUE )
  networkFlowMatrix <- dummy
  
  save( networkFlowMatrix, directCopyNumberContribution, explainedVariance, correlationStatistics, file = outputFile )
  
  if( output )
  {
      print( outputFile )
  }

}



#' Network propagation module: Compute relative patient-specific impact matrix
#'
#' This function computes a relative patient-specific impact matrix utilizing a given correlation statistics. This function requires a pre-computed basic initial relative patient-specific impact matrix computed by the function \code{\link{computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts}}. This function allows to separate the computation of the initial from the final impact matrix. See \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts}} and \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}} for more details.
#' @param patient Patient considered for impact matrix computation
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics are taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation if the impacts of the current and the previous impact matrix differ less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts}}, \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}}
#' @return Saves final impact matrix. See output for details.
#' @export
#' @examples
#'
#' #Compute impact matrix for first patient using pre-computed correlation statistics
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' computeNetworkFlowBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = 1, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
computeNetworkFlowBasedOnGivenCorStat_PatientSpecificRelativeImpacts <- function( patient, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
    computeNetworkFlow_PatientSpecificRelativeImpacts( patient = patient, dataSetName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" ), networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
}



#' Network propagation module: Compute relative patient-specific impact matrix
#'
#' This function computes a relative patient-specific impact matrix utilizing a given correlation statistics. This enables to analyze a patient (patient) in a new data set (data, dataSetName) under consideration of the predictive power of another data set (corStatDataSetName). This supports the analysis of data sets with few patients where the sample size is too small to obtain a robust estimate of the predictive power of individual genes (correlation statistics). See \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}} for more details. The correlation statistics can be pre-computed using the function \code{\link{predictGeneExpression}}.
#' @param patient Patient considered for impact matrix computation
#' @param data Data set containing gene expression and gene copy number profiles for impact computations
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics are taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param colSumsThreshold Stop iterative computation if the impacts of the current and the previous impact matrix differ less than the given threshold
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}}
#' @return Saves basic initial and resulting final impact matrices. See output for details.
#' @export
#' @examples
#'
#' #Load test and train data set and make test data set compatible to the train data set
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' gbmTFData = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#'
#' #Determine correlation statistics
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' predictGeneExpression( data = trainData, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
#' #Compute impact matrix for first patient using pre-computed correlation statistics
#' computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = 1, data = testData, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath )
#'
computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts <- function( patient, data, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, output = TRUE )
{
    ##Compute basic network flow matrix
    computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = patient, data = data, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, path = path, output = output )
    
    ##Compute final network flow matrix
    computeNetworkFlowBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = patient, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
}