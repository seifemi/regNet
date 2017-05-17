
#' Network inference module: Load network
#'
#' This function loads a specified network (networkName) and filters for non-significant regulatory links (pValCutoff).
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{learnNetwork_ParallelComputation}}
#' @return List containing the components of the network. G: network coefficient matrix (rows: response genes, columns: predictors). o: gene-specific offset. pVal: FDR-adjusted p-value matrix belonging to G. responseGenes: Response genes. responseGenesChr: Chromsomes on which the response genes are located. responseGenesLoc: Chromosomal locations of response genes. predictors: Names of preditors. 
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' loadNetworkWithFilteringForSignificantPredictors( networkName = "AS_SignatureTFs", pValCutoff = 0.01, path = projectPath )
#' 
loadNetworkWithFilteringForSignificantPredictors <- function( networkName, pValCutoff, path, output = TRUE )
{
    ##Load networks of individual runs
    networks = list()
    
    if( output )
    {
      print( "Load network" )
    }
    
    #Loads network object stored as res
    subPath = "/NetworkModel/WholeNetwork/"
    load( paste( path, subPath, networkName, "_NetworkCreator.Rout" , sep = "" ) )
    networks[[ 1 ]] = res
    
    
    ##Reduce networks based on p-value cutoff
    networkColumns = ncol( networks[[ 1 ]]$G )
    if( pValCutoff < 1 )
    {
      if( output )
      {
        print( paste( "Remove non-significant predictors with p-values >", pValCutoff ) )
      }
    
      #Determine predictors with non-significant p-value
      pValVector = as.vector( networks[[ 1 ]]$pVal )
      unstableLinks = c( which( pValVector > pValCutoff ), which( is.na( pValVector ) ) )
    
      #Network in vector representation
      dummyNetVector = as.vector( as.matrix( networks[[ 1 ]]$G ) )

      #Remove links
      if( length( unstableLinks ) > 0 )
      {
        dummyNetVector[ unstableLinks ] = 0
      }

      #Recreate network from cleaned link presentation
      networks[[ 1 ]]$G = matrix( dummyNetVector, ncol = networkColumns )
    }    
    
    return( networks[[ 1 ]] )
}
