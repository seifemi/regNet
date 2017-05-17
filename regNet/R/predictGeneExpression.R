
#' Network prediction module: Network-based prediction of gene expression data
#'
#' This function predicts expression levels of all genes in the data set (data) using the specified network (networkName). Only significant regulatory links (pValCutoff) that are not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene are included in the prediction. This function returns the correlation between the network-based predicted and originally measured expression levels of each gene across all samples in a given data set. Also the corresponding p-value specifying if a gene-specific correlation is greater than zero is returned.
#' @param data Data set containing gene expression and gene copy number profiles for prediction
#' @param dataSetName Name of the data set
#' @param networkName Name of the network
#' @param pValCutoff Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.
#' @param localGeneCutoff Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.
#' @param path Project path
#' @param output Show progress information. Default: TRUE
#' @return Saves correlation and p-value statistics. See output for details.
#' @export
#' @examples
#'
#' #Load test and train data set and make test data set compatible to the train data set
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' gbmTFData = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#' 
#' #Predict expression levels of TFs in GBMs
#' projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#' predictGeneExpression( data = gbmTFData, dataSetName = "TCGA_GBM_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )
#' 
predictGeneExpression <- function( data, dataSetName, networkName, pValCutoff, localGeneCutoff, path, output = TRUE )
{
  ##Load network and remove for non-significant predictors
  network = loadNetworkWithFilteringForSignificantPredictors( networkName = networkName, pValCutoff = pValCutoff, path = path, output = output )
  
  ##Check if genes in the data set considered for prediction are in the same order of genes as in the data set used to learn the network
  if( !identical( data$genes, network$responseGenes ) )
  {
      print( paste( "Error: Genes in the given data set '", dataSetName, "' are not in the same order as the genes in the data set used to learn the network '", networkName, "'.", sep = "" ) )
      print( "" )
      print( " - Please ensure that the given data set has exactly the same order of genes." )
      print( " - Remove genes from the given data set if they were not contained in the data set used for network learning." )
      print( " - Add missing genes with artifical measurement values of zero to the given data set." )
      print( "" )
      print( "Use the function 'makeTestDataSetCompatible' to realize this." )
      
      return()
  }
  
  ##
  ##Basic variables
  ##    
  #Number of response variables: all genes in the data set
  nrResponses = length( data$genes )

  #Number of predictor variables (gene-specific copy number + expression levels)
  nrPredictors = nrow( data$Y ) + 1

  # Number of experiments
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
  
  ##Output statistics
  corStatistics = array( NA, dim = nrResponses ) 
  corPValStatistics = array( NA, dim = nrResponses )

  
  if( output )
  {
      print( paste( "Prediction of", dataSetName, "utilizing", networkName ) )
  }
  
  for( i in 1:nrResponses )    
  {
      #print( i )
  
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
      ##Predict gene expression level of response gene
      ##
      #Predict expression level ignoring offset (later correlation evaluation is invariant to the offset)
      Yhat = specificM %*%  network$G[ i,  ]
      
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
        
        corStatistics[ i ] = dummy[[ 4 ]] #Correlation
        corPValStatistics[ i ] = dummy[[ 3 ]] #P-value
      }
  }
  
  ##Save correlation and p-value statistics
  outputTab = cbind( data$genes, corStatistics, corPValStatistics )
  colnames( outputTab ) <- c( "Gene", "Correlation", "P-Value" )
  
  subPath = "/NetworkPredictions/"
  outputFile = paste( path, subPath, dataSetName, "_PredictionOfGeneExpressionBasedOn_", networkName, "_PValueCutoff_", pValCutoff, "_LocalGeneCutoff_", localGeneCutoff ,".txt", sep = "" )
  write.table( outputTab, file = outputFile, col.names = TRUE, row.names = FALSE, quote = FALSE, dec = ".", sep = "\t" )

  if( output )
  {
      print( "Write prediction statistics:" )
      print( outputFile )
  }
  
}
