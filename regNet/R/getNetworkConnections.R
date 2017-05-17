
#' Network inference module: Get network connections
#'
#' This function determines learned links between genes in the given network (networkName). These links comprise all incoming and outgoing activator and/or inhibitor links of each individual gene. Also direct gene copy number effects are reported. Only links with p-values equal or below the given p-value cutoff (pValCutoff) are considered. For each specific gene, all incoming/outgoing links from/to genes in close chromosomal proximity up- and downstream of the considered gene are removed (localGeneCutoff).
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @return Saves a network connectivity table. See output for details.
#' @export
#' @seealso \code{\link{getNetworkConnectionsForVisualization}}
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' getNetworkConnections( networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
getNetworkConnections <- function( networkName, pValCutoff, localGeneCutoff, path, output = TRUE )
{
    ##Load network
    res = loadNetworkWithFilteringForSignificantPredictors( networkName = networkName, pValCutoff = pValCutoff, path = path, output = output )
    
    N = length( res$responseGenes )
    
    ##Remove local predictors
    if( localGeneCutoff > 0 )
    {
      if( output )
      {
        print( "Remove local predictors" )
      }
      
      for( i in 1:N )
      {
        currentResponseGene    = res$responsGenes[ i ]
        currentResponseGeneChr = res$responseGenesChr[ i ]
  
        #Determine downstream predictors
        removePredictorPos = c()
        for( v in max( 1, i - localGeneCutoff ):max( 1, i - 1 ) )
        {
          localNeighborGeneChr = res$responseGenesChr[ v ]
    
          #If local neighbor is still on same chromosome then mark it for later removement
          if( localNeighborGeneChr == currentResponseGeneChr )
          {
            removePredictorPos = c( removePredictorPos, v )
          }      
        }
  
        #Determine upstream predictors
        for( v in min( N, i + 1 ):min( N, i + localGeneCutoff ) )
        {
          localNeighborGeneChr = res$responseGenesChr[ v ]
      
          #If local neighbor is still on same chromosome then mark it for later removement
          if( localNeighborGeneChr == currentResponseGeneChr )
          {
            removePredictorPos = c( removePredictorPos, v )
          }      
        }

        res$G[ i, removePredictorPos + 1 ] = 0 ##Offset of 1 because the first predictor is the gene copy number followed by all predictor genes in the order of the response genes
      }
    }  
    
    
    ##
    ##Determine incoming and outgoing links for each gene (gene as response variable, gene as predictor variable)    
    ##
    
    if( output )
    {
      print( "Determine gene-specific network connections" )
    }
    
    connectivityTable = matrix( "", nrow = N, ncol = 11 )

    for( i in 1:N )
    {
      connectivityTable[ i, 1 ] = res$responseGenes[ i ]

      ##Incomming links (gene as response)

      #Repressive links
      negInLinks = which( res$G[ i, ] < 0 )
      negInLinkGenes = pasteIDs( ids = res$predictors[ negInLinks ], sep = ";" )

      #Activating links
      posInLinks = which( res$G[ i, ] > 0 )
      posInLinkGenes = pasteIDs( ids = res$predictors[ posInLinks ], sep = ";" )

      connectivityTable[ i, 2 ] = as.character( length( negInLinks ) )
      connectivityTable[ i, 3 ] = as.character( length( posInLinks ) )
      connectivityTable[ i, 4 ] = as.character( length( negInLinks ) + length( posInLinks ) )
      connectivityTable[ i, 8 ] = negInLinkGenes
      connectivityTable[ i, 9 ] = posInLinkGenes

      ##Outgoing links (gene as predictor)
      predictorPos = which( res$predictors == res$responseGenes[ i ] )

      #Repressive links
      negOutLinks = which( res$G[ , predictorPos ] < 0 )
      negOutLinkGenes = pasteIDs( ids = res$responseGenes[ negOutLinks ], sep = ";" )

      #Activating links
      posOutLinks = which( res$G[  , predictorPos ] > 0 )
      posOutLinkGenes = pasteIDs( ids = res$responseGenes[ posOutLinks ], sep = ";" )

      connectivityTable[ i, 5 ] = as.character( length( negOutLinks ) )
      connectivityTable[ i, 6 ] = as.character( length( posOutLinks ) )
      connectivityTable[ i, 7 ] = as.character( length( negOutLinks ) + length( posOutLinks ) )
      connectivityTable[ i, 10 ] = negOutLinkGenes
      connectivityTable[ i, 11 ] = posOutLinkGenes
    }
    
    ##Save output
    colNames = c( "Gene", "IncomingRepressorLinks", "IncomingActivatorLinks", "IncomingLinks", "OutgoingRepressorLinks", "OutgoingActivatorLinks", "OutgoingLinks", "IncomingRepressorLinkGenes", "IncomingActivatorLinkGenes",  "OutgoingRepressorLinkGenes", "OutgoingActivatorLinkGenes" )
    colnames( connectivityTable ) <- colNames

    subPath = "/NetworkConnectivity/"
    outputFile = paste( path, subPath, networkName, "_GeneSpecificNetworkConnections_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".txt" , sep = "" )
    if( output )
    {
      print( "Save network connections:" )
      print( outputFile )
    }
    write.table( connectivityTable, file = outputFile, row.names = FALSE, col.names = TRUE, quote = FALSE, dec = ".", sep = "\t" )    
}



#' Network inference module: Get network connections for visualization
#'
#' This function determines regulatory links between genes in the given network (networkName). These links comprise all incoming and outgoing activator and/or inhibitor links of each individual gene. Also direct gene copy number effects are reported. Only links with p-values equal or below the given p-value cutoff (pValCutoff) are considered. For each specific gene, all incoming/outgoing links from/to genes in close chromosomal proximity up- and downstream of the considered gene are removed (localGeneCutoff). The resulting output file can be used for network visualizations (e.g. R package igraph).
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @return Saves a network connectivity table. See output for details.
#' @export
#' @seealso \code{\link{getNetworkConnections}}
#' @examples
#'
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' getNetworkConnectionsForVisualization( networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#'
getNetworkConnectionsForVisualization <- function( networkName, pValCutoff, localGeneCutoff, path, output = TRUE )
{
    ##Load network
    res = loadNetworkWithFilteringForSignificantPredictors( networkName = networkName, pValCutoff = pValCutoff, path = path, output = output )
    
    N = length( res$responseGenes )
    
    ##Remove local predictors
    if( localGeneCutoff > 0 )
    {
      if( output )
      {
        print( "Remove local predictors" )
      }
      
      for( i in 1:N )
      {
        currentResponseGene    = res$responsGenes[ i ]
        currentResponseGeneChr = res$responseGenesChr[ i ]
  
        #Determine downstream predictors
        removePredictorPos = c()
        for( v in max( 1, i - localGeneCutoff ):max( 1, i - 1 ) )
        {
          localNeighborGeneChr = res$responseGenesChr[ v ]
      
          #If local neighbor is still on same chromosome then mark it for later removement
          if( localNeighborGeneChr == currentResponseGeneChr )
           {
             removePredictorPos = c( removePredictorPos, v )
           }      
        }
  
        #Determine upstream predictors
        for( v in min( N, i + 1 ):min( N, i + localGeneCutoff ) )
        {
          localNeighborGeneChr = res$responseGenesChr[ v ]
     
          #If local neighbor is still on same chromosome then mark it for later removement
          if( localNeighborGeneChr == currentResponseGeneChr )
          {
            removePredictorPos = c( removePredictorPos, v )
          }      
        }
      
        res$G[ i, removePredictorPos + 1 ] = 0 ##Offset of 1 because the first predictor is the gene copy number followed by all predictor genes in the order of the response genes
      }
    }  
    
    ##
    ##Determine outgoing links for each gene (gene as predictor variable)    
    ##
    
    if( output )
    {
        print( "Determine gene-specific network connections" )
    }

    
    connectivityTable = matrix( "", nrow = N, ncol = 3 )
    counter = 1
    maxCounter = N
    for( i in 1:N )
    {
        ##Outgoing links (gene as predictor)
        predictorPos = which( res$predictors == res$responseGenes[ i ] )

        #Repressive links
        negOutLinks = which( res$G[ , predictorPos ] < 0 )
        negOutLinkGenes = res$responseGenes[ negOutLinks ]
        L = length( negOutLinkGenes )
        if( L > 0 )
        {
          for( l in 1:L )
          {
            connectivityTable[ counter, 1 ] = res$responseGenes[ i ]
            connectivityTable[ counter, 2 ] = negOutLinkGenes[ l ]
            connectivityTable[ counter, 3 ] = "-1"
            counter = counter + 1

            if( counter > maxCounter )
            {
              maxCounter = maxCounter + N
              connectivityTable = rbind( connectivityTable, matrix( "", nrow = N, ncol = 3 ) )
            }

          }
        }

        #Activating links
        posOutLinks = which( res$G[  , predictorPos ] > 0 )
        posOutLinkGenes = res$responseGenes[ posOutLinks ]
        L = length( posOutLinkGenes )
        if( L > 0 )
        {
          for( l in 1:L )
          {
            connectivityTable[ counter, 1 ] = res$responseGenes[ i ]
            connectivityTable[ counter, 2 ] = posOutLinkGenes[ l ]
            connectivityTable[ counter, 3 ] = "1"
            counter = counter + 1

            if( counter > maxCounter )
            {
              maxCounter = maxCounter + N
              connectivityTable = rbind( connectivityTable, matrix( "", nrow = N, ncol = 3 ) )
            }
          }
        }
    }

    if( counter > 1 )
    {
      counter = counter - 1
      connectivityTable = connectivityTable[ 1:counter, ]
    
      ##Save output
      colNames = c( "From", "To", "Type" )
      colnames( connectivityTable ) <- colNames

      subPath = "/NetworkConnectivity/"
      outputFile = paste( path, subPath, networkName, "_GeneSpecificNetworkConnectionsForVisualization_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff, ".txt" , sep = "" )
      if( output )
      {
        print( "Save network connections:" )
        print( outputFile )
      }
      write.table( connectivityTable, file = outputFile, row.names = FALSE, col.names = TRUE, quote = FALSE, dec = ".", sep = "\t" )
   }
   else
   {
     if( output )
     {
       print( "Considered network does not contain any links." )
     }
   }
    
}



#
#Concatinates all ids entries: internal function
#
pasteIDs <- function( ids, sep )
{
    #Gene does not have a link
    N = length( ids )
    if( N == 0 ) 
    {
      return( "" )
    }
    
    #Gene does only have a single link
    if( N == 1 )
    {
      return( ids[ 1 ] )
    }
    else
    {
      #Concatinate all entries
      res = ids[ 1 ]
      for( i in 2:N )
      {
        res = paste( res, ids[ i ], sep = sep )
      }
      return( res )
    }    
}
