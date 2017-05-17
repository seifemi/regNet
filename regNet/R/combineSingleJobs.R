
#' Network inference module: Combine individual sub-networks to a global network
#'
#' This function integrates individual sub-networks into a global network.
#' @param networkName Name of network
#' @param cores Fixed number of sub-network inference problems
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @return Saves global network and overview statistics. See output for details.
#' @seealso \code{\link{learnNetwork_ParallelComputation}}, \code{\link{loadNetworkWithFilteringForSignificantPredictors}}
#' @export
#' @examples
#' 
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' combineSingleJobs( networkName = "AS_SignatureTFs", cores = 4, path = projectPath )
#'
combineSingleJobs <- function( networkName, cores, path, output = TRUE )
{
    if( output )
    {
      print( "Combine single jobs" )
    }

    G = c()
    o = c()
    pVal = c()
    responseGenes = c()
    responseGenesChr = c()
    responseGenesLoc = c()
    predictors = c()
    
    subPath = "/NetworkModel/SingleJobs/"
    cvStatistics = c()
    for( i in 1:cores )
    {
      if( output )
      {
        print( paste( i, "of", cores ) )
      }
      
      ##Load saved network object 'res'
      networkFile = paste( path, subPath, networkName, "_Cores_", cores, "_Job_", i, "_NetworkCreator.Rout",  sep = "" )
      load( file = networkFile )
      
      #Add components
      G = rbind( G, as.matrix( res$G ) )
      o = c( o, res$o )
      pVal = rbind( pVal, as.matrix( res$pVal ) )
      responseGenes = c( responseGenes, res$responseGenes )
      responseGenesChr = c( responseGenesChr, res$responseGenesChr )
      responseGenesLoc = c( responseGenesLoc, res$responseGenesLoc )
      
      if( i == 1 )
      {
        predictors   = res$predictors
      }
      
      ##Load and save cv statistics
      cvStatisticsFile = paste( path, subPath, networkName, "_Cores_", cores, "_Job_", i, "_NetworkCreator_CVStatistics.txt",  sep = "" )
      cvStat = read.delim( file = cvStatisticsFile, header = TRUE, check.names = FALSE )
      cvStatistics = rbind( cvStatistics, cvStat )      
    }
    
    #Convert network G to a sparse network
    sG <- as( G, "sparseMatrix" )
  
    #Perform FDR correction of p-values
    N = ncol( pVal )
    pValVector = as.vector( pVal )
    nonNAEntries = which( !is.na( pValVector ) )
    nonNAPValues = pValVector[ nonNAEntries ]
    fdrAdjustedPValues = p.adjust( nonNAPValues, method = "fdr" )
    pValVector[ nonNAEntries ] = fdrAdjustedPValues
    pVal = matrix( pValVector, ncol = N )  
  
    #Network G, Offset o, ...
    res = list( G = sG, o = o, pVal = pVal, responseGenes = responseGenes, responseGenesChr = responseGenesChr, responseGenesLoc = responseGenesLoc, predictors = predictors )
    
    #Save network object
    subPath = "/NetworkModel/WholeNetwork/"
    networkFile = paste( path, subPath, networkName, "_NetworkCreator.Rout",  sep = "" )
    save( res, file = networkFile )
    if( output )
    {
      print( "Save network:" )
      print( networkFile )
    }
    
    #Save cv statistics
    cvStatisticsFile = paste( path, subPath, networkName, "_NetworkCreator_CVStatistics.txt",  sep = "" )
    write.table( cvStatistics, file = cvStatisticsFile, row.names = FALSE, col.names = TRUE, sep = "\t", dec = ".", quote = FALSE )
    if( output )
    {
      print( "Save network statistics:" )
      print( cvStatisticsFile )
    }
}
