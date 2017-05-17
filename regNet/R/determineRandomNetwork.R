
#' Network inference module: Create random network instance
#'
#' This function loads a specified network (networkName), filters this network for non-significant regulatory links (pValCutoff), and finally creates a random degree-preserving instance of the network. The resulting random network is saved using a pre-defined file-naming convention including the random network instance number (randomNetworkInstance) to uniquely label the random network for later usage. The random network is compatible with all other regNet functions and can be used to generate baseline results for network predictions and network propagation.
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param randomNetworkInstance Random network instance number 
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @seealso \code{\link{loadNetworkWithFilteringForSignificantPredictors}}, \code{\link{predictGeneExpression}}, ...
#' @return Saves random network instance. See output for details.
#' @export
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' determineRandomNetworkWithFilteringForSignificantPredictors( networkName = "AS_SignatureTFs", pValCutoff = 0.01, randomNetworkInstance = 1, path = projectPath )
#' 
determineRandomNetworkWithFilteringForSignificantPredictors <- function( networkName, pValCutoff, randomNetworkInstance, path, output = TRUE )
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
    
    
    ##Reduce network based on p-value cutoff
    networkColumns = ncol( networks[[ 1 ]]$G )
    if( pValCutoff <= 1 )
    {
      if( output )
      {
        print( paste( "Remove non-significant predictors with p-values >", pValCutoff ) )
      }
    
      #Determine predictors with non-significant and significant p-value (non-significant predictors also include predictors with p-values of NA)
      pValVector = as.vector( networks[[ 1 ]]$pVal )
      unstableLinks = c( which( pValVector > pValCutoff ), which( is.na( pValVector ) ) )
      nrStableLinks   = length( which( pValVector <= pValCutoff ) )
    
      #Network in vector representation
      dummyNetVector = as.vector( as.matrix( networks[[ 1 ]]$G ) )

      #Remove links
      if( length( unstableLinks ) > 0 )
      {
        dummyNetVector[ unstableLinks ] = 0
      }

      #Recreate network from cleaned link presentation
      networks[[ 1 ]]$G = matrix( dummyNetVector, ncol = networkColumns )

      #Degree-preserving network permutation
      permutedNetwork_G = randomizeNetwork_DegreePreservingNetworkPermutation( network = networks[[ 1 ]]$G, nrLinks = nrStableLinks, output = output )

      #Create new p-value matrix: All permuted links get an artifical p-value equal to pValCutoff
      pValVector[ 1:length( pValVector ) ] = NA
      pos = which( as.vector( permutedNetwork_G ) != 0 )
      pValVector[ pos ] = pValCutoff
      pVal = matrix( pValVector, ncol = networkColumns )

      #Modify network object
      networks[[ 1 ]]$G <- as( permutedNetwork_G, "sparseMatrix" )
      networks[[ 1 ]]$pVal <- pVal
      dummyOffset = networks[[ 1 ]]$o
      dummyOffset[ 1:length( dummyOffset ) ] = 0
      networks[[ 1 ]]$o = dummyOffset
    }    
    
    
    #Save network object
    res = networks[[ 1 ]]
    networkFile = paste( path, subPath, "RandomNetwork_", randomNetworkInstance, "_PValueCutoff_", pValCutoff, "_BasedOn_", networkName, "_NetworkCreator.Rout" , sep = "" )
    save( res, file = networkFile )
    if( output )
    {
      print( networkFile )
    }
}



#
#Degree-preserving network permutation: internal function
#
#Input: network matrix, number of links
#
randomizeNetwork_DegreePreservingNetworkPermutation <- function( network, nrLinks, output = TRUE )
{
    if( output )
    {
      print( "Degree-preserving network permutation" )
    }
    
    ##Copy network object
    randomNet = network

        
    ##Get edge list
    predictorGenes = array( NA, dim = nrLinks )
    responseGenes = array( NA, dim = nrLinks )
        
    if( output )
    {
      print( "Determine trans-acting links" )
    }
        
    N = nrow( randomNet )
    counter = 1
    for( i in 1:N )
    {
      cand = setdiff( which( randomNet[ i, ] != 0 ), 1 ) ##all trans-regulators, but not local cis-CNV of a response gene
      L = length( cand )

      if( L > 0 )
      {
        for( l in 1:L )
        {
          predictorGenes[ counter ] = cand[ l ]
          responseGenes[ counter ] = i
          counter = counter + 1
        }
      }
    }
    
    cand = which( !is.na( predictorGenes ) )
    N = length( cand )
    
    if( output )
    {     
      print( paste( "   - trans-acting links:", N ) )
      print( paste( "   - cis-acting links:", nrLinks - N ) )
    }
    
    
    predictorGenes = predictorGenes[ cand ]
    responseGenes  = responseGenes[ cand ]
    
    
    ##Degree-preserving network permutation: trans-acting links
    if( output )
    {
      print( "Perform degree-preserving network permutation" )
    }
    
    sampleUniverse = 1:N
    while( length( sampleUniverse ) > 1 )
    {
      #Randomly choose two edges and exchange predictors for both response genes
      exchangeCands = sample( x = sampleUniverse, size = 2 )
      i             = exchangeCands[ 1 ]
      exchangeCand  = exchangeCands[ 2 ]

      randomNet[ responseGenes[ i ], predictorGenes[ i ] ] = 0
      randomNet[ responseGenes[ i ], predictorGenes[ exchangeCand ] ] = network[ responseGenes[ exchangeCand ], predictorGenes[ exchangeCand ] ]

      randomNet[ responseGenes[ exchangeCand ], predictorGenes[ exchangeCand ] ] = 0
      randomNet[ responseGenes[ exchangeCand ], predictorGenes[ i ] ] = network[ responseGenes[ i ], predictorGenes[ i ] ]

      #Remove re-sampled edges
      sampleUniverse = setdiff( sampleUniverse, exchangeCands )
    }
        
    ##CNV must be handled differently: Random permutation of CNV effects: cis-acting links
    randomNet[ , 1 ] = sample( randomNet[ , 1 ] )
    
    return( randomNet )        
} 
