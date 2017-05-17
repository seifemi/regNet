
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
## Compute average absolute or relative impacts
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



#' Network propagation module: Determine average absolute cohort-specific impacts
#'
#' This function determines average absolute cohort-specific impacts for each given source gene (sourceGenes) on all given target genes (targetGenes) utilizing a pre-computed cohort-specific impact matrix. The average impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the average impacts are also directly returned as matrix. Average absolute impacts allow to identify those source genes that have the greatest average impact on target genes.
#' @param sourceGenes Source genes for which average impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}, \code{\link{getAverageImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}}
#' @return Returns and saves average impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' getAverageImpacts_CohortSpecificAbsoluteImpacts( sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#' 
getAverageImpacts_CohortSpecificAbsoluteImpacts <- function( sourceGenes, targetGenes, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    if( output )
    {
      print( "Load cohort-specific network flow matrix:" )
    }
    
    ##
    ##Load network flow matrix
    ##
    ##- rows and columns contain the same number of genes in identical order
    ##- impact of a gene i on a gene j is given by the entry flowMatrix[ i, j ]
    ##
    flowMatrix = loadNetworkFlowMatrix_CohortSpecificAbsoluteImpacts( dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )

    ##Compute cohort-specific impacts
    res = getAverageImpacts_InternalFunction( flowMatrix = flowMatrix, sourceGenes = sourceGenes, targetGenes = targetGenes, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#' Network propagation module: Determine average absolute cohort-specific impacts
#'
#' This function determines average absolute cohort-specific impacts for each given source gene (sourceGenes) on all given target genes (targetGenes) utilizing a pre-computed cohort-specific impact matrix that was computed using the predictive power of a given data set (corStatDataSetName). The average impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the average impacts are also directly returned as matrix. Average absolute impacts allow to identify those source genes that have the greatest average impact on target genes.
#' @param sourceGenes Source genes for which average impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics were taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}}, \code{\link{getAverageImpacts_CohortSpecificAbsoluteImpacts}}
#' @return Returns and saves average impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' data = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#' getAverageImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getAverageImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts <- function( sourceGenes, targetGenes, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    combinedName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" )
    res = getAverageImpacts_CohortSpecificAbsoluteImpacts( sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = combinedName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#' Network propagation module: Determine average absolute patient-specific impacts
#'
#' This function determines average absolute patient-specific impacts for each given source gene (sourceGenes) on all given target genes (targetGenes) utilizing a pre-computed patient-specific impact matrix. The average impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the average impacts are also directly returned as matrix. Average absolute impacts allow to identify those source genes that have the greatest average impact on target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which average impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}, \code{\link{getAverageImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}}
#' @return Returns and saves average impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' getAverageImpacts_PatientSpecificAbsoluteImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#' 
getAverageImpacts_PatientSpecificAbsoluteImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
      if( output )
      {
        print( "Load patient-specific network flow matrix:" )
      }

      ##Load patient-specific flow matrix
      flowMatrix = loadNetworkFlowMatrix_PatientSpecificAbsoluteImpacts( patient = patient, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
      
      ##Compute patient-specific impacts
      res = getAverageImpacts_InternalFunction( flowMatrix = flowMatrix, sourceGenes = sourceGenes, targetGenes = targetGenes, path = path, outputFile = outputFile, output = output )
      
      return( res )
}



#' Network propagation module: Determine average absolute patient-specific impacts
#'
#' This function determines average absolute patient-specific impacts for each given source gene (sourceGenes) on all given target genes (targetGenes) utilizing a pre-computed patient-specific impact matrix that was computed using the predictive power of a given data set (corStatDataSetName). The average impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the average impacts are also directly returned as matrix. Average absolute impacts allow to identify those source genes that have the greatest average impact on target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which average impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics were taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}}, \code{\link{getAverageImpacts_PatientSpecificAbsoluteImpacts}}
#' @return Returns and saves average impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' data = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#' getAverageImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getAverageImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    combinedName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" )
    res = getAverageImpacts_PatientSpecificAbsoluteImpacts( patient = patient, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = combinedName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#' Network propagation module: Determine average relative patient-specific impacts
#'
#' This function determines average relative patient-specific impacts for each given source gene (sourceGenes) on all given target genes (targetGenes) utilizing a pre-computed patient-specific impact matrix. The average impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the average impacts are also directly returned as matrix. Average relative impacts account for potential activator and inhibitory contributions along the network paths from source to target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which average impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}}, \code{\link{getAverageImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts}}
#' @return Returns and saves average impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' getAverageImpacts_PatientSpecificRelativeImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getAverageImpacts_PatientSpecificRelativeImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
      if( output )
      {
        print( "Load patient-specific network flow matrix:" )
      }

      ##Load patient-specific flow matrix
      flowMatrix = loadNetworkFlowMatrix_PatientSpecificRelativeImpacts( patient = patient, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
      
      ##Compute patient-specific impacts
      res = getAverageImpacts_InternalFunction( flowMatrix = flowMatrix, sourceGenes = sourceGenes, targetGenes = targetGenes, path = path, outputFile = outputFile, output = output )
      
      return( res )
}



#' Network propagation module: Determine average relative patient-specific impacts
#'
#' This function determines average relative patient-specific impacts for each given source gene (sourceGenes) on all given target genes (targetGenes) utilizing a pre-computed patient-specific impact matrix that was computed using the predictive power of a given data set (corStatDataSetName). The average impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the average impacts are also directly returned as matrix. Average relative impacts account for potential activator and inhibitory contributions along the network paths from source to target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which average impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics were taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts}}, \code{\link{getAverageImpacts_PatientSpecificRelativeImpacts}}
#' @return Returns and saves average impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' data = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#' getAverageImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getAverageImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    combinedName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" )
    res = getAverageImpacts_PatientSpecificRelativeImpacts( patient = patient, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = combinedName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#
#General internal function for computing average absolute or relative impacts: internal function.
#
getAverageImpacts_InternalFunction <- function( flowMatrix, sourceGenes, targetGenes, path, outputFile, output )
{
    genes = colnames( flowMatrix )    
    mapPos = match( targetGenes, genes )
    
    #Restrict to those target genes that are contained in the columns of the flow matrix
    NA_dummy = which( is.na( mapPos ) )
    notNA_dummy         = which( !is.na( mapPos ) )
    mapPos              = mapPos[ notNA_dummy ]    
    nrMappedTargetGenes = length( notNA_dummy )
    
    if( output )
    {
      print( paste( "Mapped target genes:", nrMappedTargetGenes, "of", length( targetGenes ) ) )
      if( nrMappedTargetGenes != length( targetGenes ) )
      {
        print( "Not mapped target genes:" )
        print( targetGenes[ NA_dummy ] )
      }
    }
    
    #Only the columns of the target genes are relevant for determining the impact of source genes
    flowMatrix = flowMatrix[ , mapPos ]    
    mappedTargetGenes = genes[ mapPos ]

    
    if( output )
    {
      if( nrMappedTargetGenes > 1 )
      {
        dummy = colSums( flowMatrix )
      }
      else
      {
        if( nrMappedTargetGenes == 1 )
        {
          dummy = sum( flowMatrix ) 
        }
        else
        {
          if( output )
          {
            print( "No network paths from source to target genes exist!" )
          }
      
          return( NA )
        }
      }
      cand = which( dummy == 0 )      
      
      print( paste( "Mapped target genes without predictors:", length( cand ), "of", nrMappedTargetGenes ) )
      if( length( cand ) > 0 )
      {  
        print( mappedTargetGenes[ cand ] )
      }
      
      #No connections from source to target genes exist
      if( length( cand ) == nrMappedTargetGenes )
      {
        if( output )
        {
          print( "No network paths from source to target genes exist!" )
        }
        
        return( NA )
      }
    }

    
    
    ##
    ##Standardize flow matrix: Sum of direct cis and trans-effects for each survival signature gene should be one. If one then takes the average trans-impact of a mutated genes across all survival signature genes, then one obtains the average proportion.
    ##
    
    #Standardize flow matrix
    if( nrMappedTargetGenes > 1 )
    {
      for( i in 1:nrMappedTargetGenes )    
      {
        totalImpact = sum( abs( flowMatrix[ , i ] ) ) #Note: abs( .. ) is used to enable the computation of relative impacts. abs( .. ) has no influence on the computation of absolute impacts.
    
        flowMatrix[ , i ] = flowMatrix[ , i ] / totalImpact
      }
    }
    else
    {
      #nrMappedTargetGenes == 1 -> flowMatrix is a vector
      totalImpact = sum( abs( flowMatrix ) ) #See note abs( .. ) above.
      flowMatrix = flowMatrix / totalImpact
    }

    
    ##
    ##Compute the impact of source genes on target genes
    ##
    mapPos = match( sourceGenes, genes )
    
    #Restrict to those source genes that are contained in the rows of the flow matrix
    NA_dummy = which( is.na( mapPos ) )
    notNA_dummy         = which( !is.na( mapPos ) )
    mapPos              = mapPos[ notNA_dummy ]    
    mappedSourceGenes   = genes[ mapPos ]
    nrMappedSourceGenes = length( mappedSourceGenes )
    
    if( output )
    {
      print( paste( "Mapped source genes:", nrMappedSourceGenes, "of", length( sourceGenes ) ) )
      if( nrMappedSourceGenes != length( sourceGenes ) )
      {
        print( "Not mapped source genes:" )
        print( sourceGenes[ NA_dummy ] )
      }
    }
       
    
    #Only the rows of the source genes are relevant for determining the impact on target genes. The number of mapped source genes specifies with the flowMatrix is still a matrix or now a vector.
    if( nrMappedTargetGenes > 1 )
    {
      flowMatrix = flowMatrix[ mapPos, ]
    }
    else
    {
      #nrMappedTargetGenes == 1 -> flowMatrix is a vector
      flowMatrix = flowMatrix[ mapPos ]
    }
    
      
    if( nrMappedSourceGenes > 1 )
    {
      if( nrMappedTargetGenes > 1 )
      {
        averageImpacts = as.matrix( rowMeans( flowMatrix, na.rm = TRUE ) ) #na.rm = TRUE removes target genes without any predictors (all column entries of such a target gene are NaN) from the calculation
      }
      else
      {
        #Only one target gene exists
        averageImpacts = as.matrix( flowMatrix )
      }
      rownames( averageImpacts ) <- mappedSourceGenes
      colnames( averageImpacts ) <- "AverageImpactOnTargetGenes"
    }
    else
    {
      if( nrMappedSourceGenes == 1 )
      {
        #nrMappedSourceGenes == 1 -> flowMatrix is a vector
        averageImpacts = as.matrix( mean( flowMatrix, na.rm = TRUE ) ) #na.rm = TRUE removes target genes without any predictors (all column entries of such a target gene are NaN) from the calculation
        rownames( averageImpacts ) <- mappedSourceGenes
        colnames( averageImpacts ) <- "AverageImpactOnTargetGenes"
      }
      else
      {
        return( NA )
      }
    }

    
    if( output )
    {
      print( "Save output:" )
    }
    
    ##
    ##Save impacts to file
    ##
    subPath = "/NetworkPropagation/ImpactComputations/"
    outputFileWithPath = paste( path, subPath, outputFile, sep = "" )
    resTab = data.frame( mappedSourceGenes, averageImpacts, stringsAsFactors = FALSE, row.names = NULL )
    colnames( resTab ) <- c( "SourceGene", "AverageImpactOnTargetGenes" )
    write.table( resTab, file = outputFileWithPath, col.names = TRUE, row.names = FALSE, quote = FALSE, dec = ".", sep = "\t" )
    
  
    if( output )
    {
      print( outputFileWithPath )
    }
    
    
    ##
    ##Return computed impacts
    ##
    return( resTab )
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
## Compute absolute or relative impacts
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



#' Network propagation module: Determine absolute cohort-specific impacts
#'
#' This function determines absolute cohort-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed cohort-specific impact matrix. The impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the impacts are also directly returned as matrix. Absolute impacts allow to identify those source genes that have the greatest impact on target genes.
#' @param sourceGenes Source genes for which impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}, \code{\link{getImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}}
#' @return Returns and saves impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' getImpacts_CohortSpecificAbsoluteImpacts( sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#' 
getImpacts_CohortSpecificAbsoluteImpacts <- function( sourceGenes, targetGenes, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    if( output )
    {
      print( "Load cohort-specific network flow matrix:" )
    }
    
    ##
    ##Load network flow matrix
    ##
    ##- rows and columns contain the same number of genes in identical order
    ##- impact of a gene i on a gene j is given by the entry flowMatrix[ i, j ]
    ##
    flowMatrix = loadNetworkFlowMatrix_CohortSpecificAbsoluteImpacts( dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )

    ##Compute cohort-specific impacts
    res = getImpacts_InternalFunction( flowMatrix = flowMatrix, sourceGenes = sourceGenes, targetGenes = targetGenes, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#' Network propagation module: Determine absolute cohort-specific impacts
#'
#' This function determines absolute cohort-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed cohort-specific impact matrix that was computed using the predictive power of a given data set (corStatDataSetName). The impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the impacts are also directly returned as matrix. Absolute impacts allow to identify those source genes that have the greatest impact on target genes.
#' @param sourceGenes Source genes for which impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics were taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}}, \code{\link{getImpacts_CohortSpecificAbsoluteImpacts}}
#' @return Returns and saves impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' data = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#' getImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts <- function( sourceGenes, targetGenes, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    combinedName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" )
    res = getImpacts_CohortSpecificAbsoluteImpacts( sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = combinedName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#' Network propagation module: Determine absolute patient-specific impacts
#'
#' This function determines absolute patient-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed patient-specific impact matrix. The impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the impacts are also directly returned as matrix. Absolute impacts allow to identify those source genes that have the greatest impact on target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}, \code{\link{getImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}}
#' @return Returns and saves impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' getImpacts_PatientSpecificAbsoluteImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#' 
getImpacts_PatientSpecificAbsoluteImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
      if( output )
      {
        print( "Load patient-specific network flow matrix:" )
      }

      ##Load patient-specific flow matrix
      flowMatrix = loadNetworkFlowMatrix_PatientSpecificAbsoluteImpacts( patient = patient, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
      
      ##Compute patient-specific impacts
      res = getImpacts_InternalFunction( flowMatrix = flowMatrix, sourceGenes = sourceGenes, targetGenes = targetGenes, path = path, outputFile = outputFile, output = output )
      
      return( res )
}



#' Network propagation module: Determine absolute patient-specific impacts
#'
#' This function determines absolute patient-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed patient-specific impact matrix that was computed using the predictive power of a given data set (corStatDataSetName). The impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the impacts are also directly returned as matrix. Absolute impacts allow to identify those source genes that have the greatest impact on target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics were taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}}, \code{\link{getImpacts_PatientSpecificAbsoluteImpacts}}
#' @return Returns and saves impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' data = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#' getImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    combinedName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" )
    res = getImpacts_PatientSpecificAbsoluteImpacts( patient = patient, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = combinedName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#' Network propagation module: Determine relative patient-specific impacts
#'
#' This function determines relative patient-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed patient-specific impact matrix. The impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the impacts are also directly returned as matrix. Relative impacts account for potential activator and inhibitory contributions along the network paths from source to target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}}, \code{\link{getImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts}}
#' @return Returns and saves impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' getImpacts_PatientSpecificRelativeImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getImpacts_PatientSpecificRelativeImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
      if( output )
      {
        print( "Load patient-specific network flow matrix:" )
      }

      ##Load patient-specific flow matrix
      flowMatrix = loadNetworkFlowMatrix_PatientSpecificRelativeImpacts( patient = patient, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
      
      ##Compute patient-specific impacts
      res = getImpacts_InternalFunction( flowMatrix = flowMatrix, sourceGenes = sourceGenes, targetGenes = targetGenes, path = path, outputFile = outputFile, output = output )
      
      return( res )
}



#' Network propagation module: Determine relative patient-specific impacts
#'
#' This function determines relative patient-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed patient-specific impact matrix that was computed using the predictive power of a given data set (corStatDataSetName). The impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the impacts are also directly returned as matrix. Relative impacts account for potential activator and inhibitory contributions along the network paths from source to target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics were taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts}}, \code{\link{getImpacts_PatientSpecificRelativeImpacts}}
#' @return Returns and saves impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' data = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#' getImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    combinedName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" )
    res = getImpacts_PatientSpecificRelativeImpacts( patient = patient, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = combinedName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#
#General internal function for computing average absolute or relative impacts
#
getImpacts_InternalFunction <- function( flowMatrix, sourceGenes, targetGenes, path, outputFile, output )
{
    genes = colnames( flowMatrix )    
    mapPos = match( targetGenes, genes )
    
    #Restrict to those target genes that are contained in the columns of the flow matrix
    NA_dummy = which( is.na( mapPos ) )
    notNA_dummy         = which( !is.na( mapPos ) )
    mapPos              = mapPos[ notNA_dummy ]    
    nrMappedTargetGenes = length( notNA_dummy )
    
    if( output )
    {
      print( paste( "Mapped target genes:", nrMappedTargetGenes, "of", length( targetGenes ) ) )
      if( nrMappedTargetGenes != length( targetGenes ) )
      {
        print( "Not mapped target genes:" )
        print( targetGenes[ NA_dummy ] )
      }
    }
    
    #Only the columns of the target genes are relevant for determining the impact of source genes
    flowMatrix = flowMatrix[ , mapPos ]    
    mappedTargetGenes = genes[ mapPos ]

    
    if( output )
    {
      if( nrMappedTargetGenes > 1 )
      {
        dummy = colSums( flowMatrix )
      }
      else
      {
        if( nrMappedTargetGenes == 1 )
        {
          dummy = sum( flowMatrix ) 
        }
        else
        {
          if( output )
          {
            print( "No network paths from source to target genes exist!" )
          }
    
          return( NA )
        }
      }
      cand = which( dummy == 0 )      
      
      print( paste( "Mapped target genes without predictors:", length( cand ), "of", nrMappedTargetGenes ) )
      if( length( cand ) > 0 )
      {  
        print( mappedTargetGenes[ cand ] )
      }
      
      #No connections from source to target genes exist
      if( length( cand ) == nrMappedTargetGenes )
      {
        if( output )
        {
          print( "No network paths from source to target genes exist!" )
        }
        
        return( NA )
      }
    }

    
    
    ##
    ##Standardize flow matrix: Sum of direct cis and trans-effects for each survival signature gene should be one. If one then takes the average trans-impact of a mutated genes across all survival signature genes, then one obtains the average proportion.
    ##
    
    #Standardize flow matrix
    if( nrMappedTargetGenes > 1 )
    {
      for( i in 1:nrMappedTargetGenes )    
      {
        totalImpact = sum( abs( flowMatrix[ , i ] ) ) #Note: abs( .. ) is used to enable the computation of relative impacts. abs( .. ) has no influence on the computation of absolute impacts.
    
        flowMatrix[ , i ] = flowMatrix[ , i ] / totalImpact
      }
    }
    else
    {
      #nrMappedTargetGenes == 1 -> flowMatrix is a vector      
      totalImpact = sum( abs( flowMatrix ) ) #See note abs( .. ) above.
      flowMatrix = flowMatrix / totalImpact
    }

    
    ##
    ##Get impacts of source genes on target genes
    ##
    mapPos = match( sourceGenes, genes )
    
    #Restrict to those source genes that are contained in the rows of the flow matrix
    NA_dummy = which( is.na( mapPos ) )
    notNA_dummy         = which( !is.na( mapPos ) )
    mapPos              = mapPos[ notNA_dummy ]    
    mappedSourceGenes   = genes[ mapPos ]
    nrMappedSourceGenes = length( mappedSourceGenes )
    
    if( output )
    {
      print( paste( "Mapped source genes:", nrMappedSourceGenes, "of", length( sourceGenes ) ) )
      if( nrMappedSourceGenes != length( sourceGenes ) )
      {
        print( "Not mapped source genes:" )
        print( sourceGenes[ NA_dummy ] )
      }
    }
    
    if( nrMappedSourceGenes == 0 )
    {
      if( output )
      {
        print( "No network paths from source to target genes exist!" )
      }

      return( NA )
    }
    
    #Only the rows of the source genes are relevant for determining the impact on target genes. The number of mapped source genes specifies with the flowMatrix is still a matrix or now a vector.
    if( nrMappedTargetGenes > 1 )
    {
      flowMatrix = flowMatrix[ mapPos, ]
    }
    else
    {
      #nrMappedTargetGenes == 1 -> flowMatrix is a vector
      flowMatrix = flowMatrix[ mapPos ]
    }
    
       
    if( output )
    {
      print( "Save output:" )
    }
    
    ##
    ##Save impacts to file
    ##
    subPath = "/NetworkPropagation/ImpactComputations/"
    outputFileWithPath = paste( path, subPath, outputFile, sep = "" )
    resTab = data.frame( mappedSourceGenes, flowMatrix, stringsAsFactors = FALSE, row.names = NULL )
    colnames( resTab ) <- c( "SourceGene->TargetGene", mappedTargetGenes )
    write.table( resTab, file = outputFileWithPath, col.names = TRUE, row.names = FALSE, quote = FALSE, dec = ".", sep = "\t" )
    
  
    if( output )
    {
      print( outputFileWithPath )
    }
    
    ##
    ##Return computed impacts
    ##
    
    return( resTab )
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
## Get raw absolute or relative impacts
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



#' Network propagation module: Determine raw absolute cohort-specific impacts
#'
#' This function determines raw absolute cohort-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed cohort-specific impact matrix. Raw impacts represent the entries of the impact matrix without furhter standardizing each column sum of the impact matrix to one. The raw impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the raw impacts are also directly returned as matrix. Raw absolute impacts allow to identify those source genes that have the greatest impact on target genes.
#' @param sourceGenes Source genes for which raw impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}, \code{\link{getRawImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}}
#' @return Returns and saves raw impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' getRawImpacts_CohortSpecificAbsoluteImpacts( sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#' 
getRawImpacts_CohortSpecificAbsoluteImpacts <- function( sourceGenes, targetGenes, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    if( output )
    {
      print( "Load cohort-specific network flow matrix:" )
    }
    
    ##
    ##Load network flow matrix
    ##
    ##- rows and columns contain the same number of genes in identical order
    ##- impact of a gene i on a gene j is given by the entry flowMatrix[ i, j ]
    ##
    flowMatrix = loadNetworkFlowMatrix_CohortSpecificAbsoluteImpacts( dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )

    ##Compute cohort-specific impacts
    res = getRawImpacts_InternalFunction( flowMatrix = flowMatrix, sourceGenes = sourceGenes, targetGenes = targetGenes, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#' Network propagation module: Determine raw absolute cohort-specific impacts
#'
#' This function determines raw absolute cohort-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed cohort-specific impact matrix that was computed using the predictive power of a given data set (corStatDataSetName). Raw impacts represent the entries of the impact matrix without furhter standardizing each column sum of the impact matrix to one. The raw impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the raw impacts are also directly returned as matrix. Raw absolute impacts allow to identify those source genes that have the greatest raw impact on target genes.
#' @param sourceGenes Source genes for which raw impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics were taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts}}, \code{\link{getRawImpacts_CohortSpecificAbsoluteImpacts}}
#' @return Returns and saves raw impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' data = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#' getRawImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getRawImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts <- function( sourceGenes, targetGenes, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    combinedName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" )
    res = getRawImpacts_CohortSpecificAbsoluteImpacts( sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = combinedName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#' Network propagation module: Determine raw absolute patient-specific impacts
#'
#' This function determines raw absolute patient-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed patient-specific impact matrix. Raw impacts represent the entries of the impact matrix without furhter standardizing each column sum of the impact matrix to one. The raw impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the raw impacts are also directly returned as matrix. Raw absolute impacts allow to identify those source genes that have the greatest raw impact on target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which raw impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts}}, \code{\link{getRawImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}}
#' @return Returns and saves raw impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' getRawImpacts_PatientSpecificAbsoluteImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#' 
getRawImpacts_PatientSpecificAbsoluteImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
      if( output )
      {
        print( "Load patient-specific network flow matrix:" )
      }

      ##Load patient-specific flow matrix
      flowMatrix = loadNetworkFlowMatrix_PatientSpecificAbsoluteImpacts( patient = patient, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
      
      ##Compute patient-specific impacts
      res = getRawImpacts_InternalFunction( flowMatrix = flowMatrix, sourceGenes = sourceGenes, targetGenes = targetGenes, path = path, outputFile = outputFile, output = output )
      
      return( res )
}



#' Network propagation module: Determine raw absolute patient-specific impacts
#'
#' This function determines raw absolute patient-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed patient-specific impact matrix that was computed using the predictive power of a given data set (corStatDataSetName). Raw impacts represent the entries of the impact matrix without furhter standardizing each column sum of the impact matrix to one. The raw impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the raw impacts are also directly returned as matrix. Raw absolute impacts allow to identify those source genes that have the greatest raw impact on target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which raw impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics were taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeBasicNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts}}, \code{\link{getRawImpacts_PatientSpecificAbsoluteImpacts}}
#' @return Returns and saves raw impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' data = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#' getRawImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getRawImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    combinedName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" )
    res = getRawImpacts_PatientSpecificAbsoluteImpacts( patient = patient, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = combinedName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#' Network propagation module: Determine raw relative patient-specific impacts
#'
#' This function determines raw relative patient-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed patient-specific impact matrix. Raw impacts represent the entries of the impact matrix without furhter standardizing each absolute column sum of the impact matrix to one. The raw impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the raw impacts are also directly returned as matrix. Raw relative impacts account for potential activator and inhibitory contributions along the network paths from source to target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which raw impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrix_PatientSpecificRelativeImpacts}}, \code{\link{getRawImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts}}
#' @return Returns and saves raw impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' getRawImpacts_PatientSpecificRelativeImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getRawImpacts_PatientSpecificRelativeImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
      if( output )
      {
        print( "Load patient-specific network flow matrix:" )
      }

      ##Load patient-specific flow matrix
      flowMatrix = loadNetworkFlowMatrix_PatientSpecificRelativeImpacts( patient = patient, dataSetName = dataSetName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, output = output )
      
      ##Compute patient-specific impacts
      res = getRawImpacts_InternalFunction( flowMatrix = flowMatrix, sourceGenes = sourceGenes, targetGenes = targetGenes, path = path, outputFile = outputFile, output = output )
      
      return( res )
}



#' Network propagation module: Determine raw relative patient-specific impacts
#'
#' This function determines raw relative patient-specific impacts for each given source gene (sourceGenes) on each given target gene (targetGenes) utilizing a pre-computed patient-specific impact matrix that was computed using the predictive power of a given data set (corStatDataSetName). Raw impacts represent the entries of the impact matrix without furhter standardizing each absolute column sum of the impact matrix to one. The raw impacts are saved to a user-defined file (outputFile) in the pre-defined regNet folder structure. In addition, the raw impacts are also directly returned as matrix. Raw relative impacts account for potential activator and inhibitory contributions along the network paths from source to target genes.
#' @param patient Patient for which impacts are computed
#' @param sourceGenes Source genes for which raw impacts on target genes are determined
#' @param targetGenes Target genes
#' @param dataSetName Name of the data set
#' @param corStatDataSetName Name of the data set from which pre-computed correlation statistics were taken
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links used for pre-computed impact matrix
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity used for pre-computed impact matrix
#' @param colSumsThreshold Stop criterion used for pre-computed impact matrix
#' @param path Project path
#' @param outputFile Name of output file
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts}}, \code{\link{getRawImpacts_PatientSpecificRelativeImpacts}}
#' @return Returns and saves raw impacts. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' data = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#' getRawImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = 1, sourceGenes = c( "TBR1", "CCNL2" ), targetGenes = data$genes, dataSetName = "TCGA_GBM_SignatureTFs", corStatDataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = "AvgImpactsOfSelectedRegulators.txt" )
#'
getRawImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts <- function( patient, sourceGenes, targetGenes, dataSetName, corStatDataSetName, networkName, pValCutoff, localGeneCutoff, colSumsThreshold, path, outputFile, output = TRUE )
{
    combinedName = paste( dataSetName, "_UsingCorStat_", corStatDataSetName, sep = "" )
    res = getRawImpacts_PatientSpecificRelativeImpacts( patient = patient, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = combinedName, networkName = networkName, pValCutoff = pValCutoff, localGeneCutoff = localGeneCutoff, colSumsThreshold = colSumsThreshold, path = path, outputFile = outputFile, output = output )
    
    return( res )
}



#
#General internal function for computing average absolute or relative impacts
#
getRawImpacts_InternalFunction <- function( flowMatrix, sourceGenes, targetGenes, path, outputFile, output )
{
    genes = colnames( flowMatrix )    
    mapPos = match( targetGenes, genes )
    
    #Restrict to those target genes that are contained in the columns of the flow matrix
    NA_dummy = which( is.na( mapPos ) )
    notNA_dummy         = which( !is.na( mapPos ) )
    mapPos              = mapPos[ notNA_dummy ]    
    nrMappedTargetGenes = length( notNA_dummy )
    
    if( output )
    {
      print( paste( "Mapped target genes:", nrMappedTargetGenes, "of", length( targetGenes ) ) )
      if( nrMappedTargetGenes != length( targetGenes ) )
      {
        print( "Not mapped target genes:" )
        print( targetGenes[ NA_dummy ] )
      }
    }
    
    #Only the columns of the target genes are relevant for determining the impact of source genes
    flowMatrix = flowMatrix[ , mapPos ]    
    mappedTargetGenes = genes[ mapPos ]

    
    if( output )
    {
      if( nrMappedTargetGenes > 1 )
      {
        dummy = colSums( flowMatrix )
      }
      else
      {
        if( nrMappedTargetGenes == 1 )
        {
          dummy = sum( flowMatrix ) 
        }
        else
        {
          if( output )
          {
            print( "No network paths from source to target genes exist!" )
          }
  
          return( NA )
        }
      }
      cand = which( dummy == 0 )      
      
      print( paste( "Mapped target genes without predictors:", length( cand ), "of", nrMappedTargetGenes ) )
      if( length( cand ) > 0 )
      {  
        print( mappedTargetGenes[ cand ] )
      }
      
      #No connections from source to target genes exist
      if( length( cand ) == nrMappedTargetGenes )
      {
        if( output )
        {
          print( "No network paths from source to target genes exist!" )
        }
        
        return( NA )
      }
    }

       
        
    ##
    ##Get impacts of source genes on target genes
    ##
    mapPos = match( sourceGenes, genes )
    
    #Restrict to those source genes that are contained in the rows of the flow matrix
    NA_dummy = which( is.na( mapPos ) )
    notNA_dummy         = which( !is.na( mapPos ) )
    mapPos              = mapPos[ notNA_dummy ]    
    mappedSourceGenes   = genes[ mapPos ]
    nrMappedSourceGenes = length( mappedSourceGenes )
    
    if( output )
    {
      print( paste( "Mapped source genes:", nrMappedSourceGenes, "of", length( sourceGenes ) ) )
      if( nrMappedSourceGenes != length( sourceGenes ) )
      {
        print( "Not mapped source genes:" )
        print( sourceGenes[ NA_dummy ] )
      }
    }
    
    if( nrMappedSourceGenes == 0 )
    {
      if( output )
      {
        print( "No network paths from source to target genes exist!" )
      }

      return( NA )
    }
    
    #Only the rows of the source genes are relevant for determining the impact on target genes. The number of mapped source genes specifies with the flowMatrix is still a matrix or now a vector.
    if( nrMappedTargetGenes > 1 )
    {
      flowMatrix = flowMatrix[ mapPos, ]
    }
    else
    {
      #nrMappedTargetGenes == 1 -> flowMatrix is a vector
      flowMatrix = flowMatrix[ mapPos ]
    }
    
       
    if( output )
    {
      print( "Save output:" )
    }
    
    ##
    ##Save impacts to file
    ##
    subPath = "/NetworkPropagation/ImpactComputations/"
    outputFileWithPath = paste( path, subPath, outputFile, sep = "" )
    resTab = data.frame( mappedSourceGenes, flowMatrix, stringsAsFactors = FALSE, row.names = NULL )
    colnames( resTab ) <- c( "SourceGene->TargetGene", mappedTargetGenes )
    write.table( resTab, file = outputFileWithPath, col.names = TRUE, row.names = FALSE, quote = FALSE, dec = ".", sep = "\t" )
    
  
    if( output )
    {
      print( outputFileWithPath )
    }
    
    
    ##
    ##Return computed impacts
    ##
    
    return( resTab )
}
