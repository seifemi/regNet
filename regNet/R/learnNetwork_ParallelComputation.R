
#' Network inference module: Define sub-network inference problems
#'
#' This is an internal function that divides the total number of genes (nrResponses) in a data set into a fixed number (cores) of equal-sized intervals of genes. This is required to separate the global network inference problem into independent sub-network inference problems.
#' @param cores Fixed number of sub-network inference problems
#' @param nrResponses Total number of genes in a data set
#' @return Array of (nearly) equal-sized intervals of genes
#' @export
#' @examples
#' 
#' range = getComputationRange( cores = 400, nrResponses = 15811 )
#'
getComputationRange <- function( cores, nrResponses )
{
    return( round( seq( 0, nrResponses, length.out = cores + 1 ), dig = 0 ) )
}



#' Network inference module: Learn regulatory network
#'
#' This function provides the basics for the inference of a regulatory network from the given data set (data). Network inference is usually very time consuming and therefore separated into independent sub-network inference problems (cores). Thus, this function learns putative regulators for all genes covered by the specific sub-network inference problem (job). Genes coverd by the job-specific sub-network inference problem are given by range[ job ] + 1 ... range[ job + 1 ] specified by \code{\link{getComputationRange}}. All individual job-specific sub-network inferrence problems (job from 1 to cores) can be solved in parallel or step-by-step depending on the available hardware. The global regulatory network can be constructed using \code{\link{combineSingleJobs}} after all sub-network inference problems have been solved. 
#' @param data Data set containing gene expression and gene copy number profiles for network inference
#' @param networkName Name of the network
#' @param cores Fixed number of sub-network inference problems
#' @param job Sub-network inference problem to be solved
#' @param path Project path
#' @param nfolds Number of sub-samples considered in the cross-validation step used to determine an optimal lambda for a specific gene. Default: 10
#' @param cvReplicates Number of repeats of nfold-cross-validations used to determine an optimal lambda for a specific gene. Default: 10
#' @param output Show progress information. Default: TRUE
#' @return Saves sub-network and cross-validation statistics files. Adds runtime entry to global runtime table. See output for details.
#' @seealso \code{\link{getComputationRange}}, \code{\link{combineSingleJobs}}, \code{\link{loadNetworkWithFilteringForSignificantPredictors}}
#' @export
#' @examples
#' 
#' #
#' #Solve the sub-network inference problems for a given data set. The following toy example does the network inference sequentially within some minutes. For large data sets, network inference should be done on a compute server with parallel execution of individual \code{learnNetwork_ParallelComputation} calls.
#' #
#'
#' #Load data set
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#'
#' #Creater basic folder structure
#' projectPath = createBasicFolderStructure( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#'
#' #Solve individual sub-network inference problems
#' learnNetwork_ParallelComputation( data = data, networkName = "AS_SignatureTFs", cores = 4, job = 1, path = projectPath )
#' learnNetwork_ParallelComputation( data = data, networkName = "AS_SignatureTFs", cores = 4, job = 2, path = projectPath )
#' learnNetwork_ParallelComputation( data = data, networkName = "AS_SignatureTFs", cores = 4, job = 3, path = projectPath )
#' learnNetwork_ParallelComputation( data = data, networkName = "AS_SignatureTFs", cores = 4, job = 4, path = projectPath )
#'
#' #Combine individual sub-networks to a global network
#' combineSingleJobs( networkName = "AS_SignatureTFs", cores = 4, path = projectPath )
#'
learnNetwork_ParallelComputation <- function( data, networkName, cores, job, path, nfolds = 10, cvReplicates = 10, output = TRUE )
{
    #Start date and time
    startTime = Sys.time()
    if( output )
    {
	print( startTime )
    }

   
   
    #alpha = 1: lasso
    #alpha = 0: ridge
    #0 < alpha < 1: elastic net
    alpha = 1

    ##
    ##Determine the computation range based on the number of available cores and the specified job
    ##
    responseRanges = getComputationRange( cores = cores, nrResponses = nrow( data$Y ) )
    selectedResponseRange = ( responseRanges[ job ] + 1 ):responseRanges[ job + 1 ]

    ##
    ##Basic variables
    ##    
    #Number of response variables
    nrResponses = length( selectedResponseRange )


    #Number of predictor variables (gene-specific copy number + expression levels)
    nrPredictors = nrow( data$Y ) + 1

    #Number of experiments
    T = ncol( data$U )
   
   
    ##
    ##Basic matrices: Y and U only contain data of considered response genes; M must be complete
    ##   
    #Response data
    #Gene expression levels of response genes: Y matrix -> number of experiments x number of response genes
    Y   = t( data.matrix( data$Y[ selectedResponseRange, ] ) )
    
    #Predictor data
    #Copy numbers of response genes: U matrix -> number of experiments x number of response genes
    U   = t( data.matrix( data$U[ selectedResponseRange, ] ) )      
    #Gene expression levels of predictor genes: M matrix -> number of experiments x number of predictor genes
    M   = t( data.matrix( data$Y ) )
    
    
    ##
    ##Network basics
    ##

    #
    #Network basics (G: interaction matrix, o: offset vector, pVal: p-value matrix)
    #
    G    = matrix( 0, nrow = nrResponses, ncol = nrPredictors )    
    o    = vector( mode = "numeric", length = nrResponses )
    pVal = matrix( NA, nrow = nrResponses, ncol = nrPredictors )
    
    #
    #Quality of fit
    #
    #Mean cross-validation error    
    mcv = vector( mode = "numeric", length = nrResponses )
    #Standard deviation of mean cross-validation error
    sdmcv = vector( mode = "numeric", length = nrResponses )
    #Number of non-zero parameters
    nzero = vector( mode = "numeric", length = nrResponses )
    #Optimal lambda minimizing mean cross-validation error
    optLambda = vector( mode = "numeric", length = nrResponses )
    
    #
    #CV Statistics
    #
    cvStatistics = matrix( 0, nrow = nrResponses, ncol = 7 )
    colnames( cvStatistics ) <- c( "MeanOptimalLambda", "SdOptimalLambda", "MeanOfMeanCVError", "SdOfMeanCVError", "MeanOptimalLambdaModel_ResidualSumOfSquares", "MeanOptimalLambdaModel_RootMeanSquareError", "SelectedPredictors" )
    
    for( i in 1:nrResponses )    
    {
      if( output )
      {
        print( paste( data$genes[ selectedResponseRange[ i ] ], "---", i, "of", nrResponses ) )
      }

      ##
      ##Response vector
      ##
      y = Y[ , i ]

      ##
      ##Predictor matrix: response gene specific copy number + predictors
      ##
      specificM = M	

      #Remove currently considered response gene from predictor matrix
      specificM[ , selectedResponseRange[ i ] ] = 0

      #Add copy number of response gene as a predictor
      specificM = cbind( U[ , i ], specificM )

      ##
      ##Determine optimal mean lambda by several CV runs
      ##
      cvRep_mcv   = vector( mode = "numeric", length = cvReplicates )
      cvRep_sdmcv = vector( mode = "numeric", length = cvReplicates )
      cvRep_nzero = vector( mode = "numeric", length = cvReplicates )
      cvRep_optLambda = vector( mode = "numeric", length = cvReplicates )

      for( r in 1:cvReplicates )
      {    
        #Compute coefficients
        cvRes = cv.glmnet( x = specificM, y = y, family = "gaussian", alpha = alpha, nfolds = nfolds, type.measure = "mse" )
    
        cvRep_optLambda[ r ] = cvRes$lambda.min
        pos = which( cvRes$lambda == cvRep_optLambda[ r ] )
        cvRep_mcv[ r ]   = cvRes$cvm[ pos ]
        cvRep_sdmcv[ r ] = cvRes$cvsd[ pos ]
        cvRep_nzero[ r ] = cvRes$nzero[ pos ]
      }

      ##
      ##Determine coefficients using the average optimal lambda
      ##
      optLambda[ i ] = mean( cvRep_optLambda )
      mcv[ i ]       = mean( cvRep_mcv )
      sdmcv[ i ]     = sd( cvRep_mcv ) #Standard deviation of mcv (different from mean of cvRep_sdmcv representing differences for each previous individual cv run)

      res = glmnet( x = specificM, y = y, family = "gaussian", alpha = alpha, lambda = optLambda[ i ] )
      nzero[ i ] = res$df[ 1 ]

      if( output )
      {
        print( paste( "Optimal Lambda:", optLambda[ i ] ) )
        print( paste( "Sd Optimal Lambda:", sd( cvRep_optLambda ) ) )
        print( paste( "Selected Predictors:", nzero[ i ] ) )
      }

      #Estimate Yhat based on coefficients and determine error
      Yhat = specificM %*% as.vector( res$beta[ , 1 ] ) + res$a0[ 1 ]
      rss  = sum( ( Y[ , i ] - Yhat )^2 )
      rmse = ( rss / T )^0.5

      if( output )
      {
        print( paste( "RSS:", rss ) )
        print( paste( "RMSE:", rmse ) )
      }

      #Set coefficients for gene i and set intercept
      G[ i, ]= as.vector( res$beta[ , 1 ] )
      o[ i ] = res$a0[ 1 ]

      #
      #Add to cvStatistics: c( "MeanOptimalLambda", "SdOptimalLambda", "MeanOfMeanCVError", "SdOfMeanCVError", "MeanOptimalLambdaModel_ResidualSumOfSquares", "MeanOptimalLambdaModel_RootMeanSquareError", "SelectedPredictors" )
      #
      cvStatistics[ i, 1 ] = optLambda[ i ]
      cvStatistics[ i, 2 ] = sd( cvRep_optLambda )
      cvStatistics[ i, 3 ] = mcv[ i ]
      cvStatistics[ i, 4 ] = sdmcv[ i ]
      cvStatistics[ i, 5 ] = rss
      cvStatistics[ i, 6 ] = rmse
      cvStatistics[ i, 7 ] = nzero[ i ]


      ##
      ##Determine p-values of selected predictors
      ##
      if( nzero[ i ] > 0 )
      {
        selectedPredictors = which( G[ i, ] != 0 )
        specificM_reducedToSelectedPredictors = as.matrix( specificM[ , selectedPredictors ] ) #as.matrix: necessary otherwise an error occurs if only one predictor has been selected
 
        #Compute lasso solution paths from zero to least squares fit
        res    = tryCatch( lars( x = specificM_reducedToSelectedPredictors, y = y, type = "lasso" ), error = function( e ){ print( e ); return( lars( x = specificM_reducedToSelectedPredictors, y = y, type = "lasso", max.steps = 20 ) ) } )  
  
        #Determine siginficance of predictors (predictor number, drop in covariance, p-value)
        resSig = tryCatch( as.matrix( covTest( res, x = specificM_reducedToSelectedPredictors, y = y )$results ), error = function( e ){ print( e ); print( paste( "covTestError:", data$responseGenes[ selectedResponseRange[ i ] ] ) ); return( NULL ) } )
  
        #Catch covTest errros: colinear matrix x or problems in predict.lars
        if( !is.null( resSig ) )
        {
          #Remove entries of predictors that where droped from the model after entering (note: if a predictor is entered more than once only its last p-value will be entered into pVal matrix
          if( nzero[ i ] == 1 )
          {
            resSig = t( as.matrix( resSig[ which( !is.na( resSig[ , 3 ] ) ), ] ) )
          }
          else
          {
            resSig = resSig[ which( !is.na( resSig[ , 3 ] ) ), ]
          }
  
          #Matrix structure could be lost after NA removement: fixed to have again a matrix
          if( is.null( nrow( resSig ) ) )
          {
            resSig = t( as.matrix( resSig ) )
          }
  
          #Set p-values for selected predictors (note: not all predictors must have been selected by covTest, some of them may not have been used and will still have NA in pVal
          pVal[ i, selectedPredictors[ resSig[ , 1 ] ] ] = resSig[ , 3 ]
        }
  
        if( output )
        {
          print( "Predictor p-value computation done" )
        }
  
      }

    }
    
    
    #Convert network G to a sparse network
    sG <- as( G, "sparseMatrix" )
    
    #Network G, Offset o, ...
    res = list( G = sG, o = o, pVal = pVal, responseGenes = data$genes[ selectedResponseRange ], responseGenesChr = data$chr[ selectedResponseRange ], responseGenesLoc = data$loc[ selectedResponseRange ], predictors = c( "GeneSpecificCopyNumber", data$genes )  )
    
    #Save network object
    subPath = "/NetworkModel/SingleJobs/"
    networkFile = paste( path, subPath, networkName, "_Cores_", cores, "_Job_", job, "_NetworkCreator.Rout",  sep = "" )
    save( res, file = networkFile )
    if( output )
    {
      print( networkFile )
    }


    #Save cv-statistics
    cvStatisticsFile = paste( path, subPath, networkName, "_Cores_", cores, "_Job_", job, "_NetworkCreator_CVStatistics.txt",  sep = "" )
    if( output )
    {
      print( cvStatisticsFile )
    }
    Genes = data$genes[ selectedResponseRange[ 1:nrResponses ] ]
    write.table( cbind( Genes, cvStatistics ), file = cvStatisticsFile, row.names = FALSE, col.names = TRUE, sep = "\t", dec = ".", quote = FALSE )
        
    #Stop date and time
    stopTime = Sys.time()
    if( output )
    {
      print( stopTime )
    }
    
    #Complete runtime
    runTime  = as.numeric( difftime( stopTime, startTime, units = "secs" ) )
    if( output )
    {
      print( paste( "Runtime:", runTime, "sec" ) )
    }

    #Add runtime of job to the corresponding output file
    outTable = cbind( job, runTime )
    colnames( outTable ) <- c( "Job", "RuntimeInSecs" )
    subPath = "/NetworkModel/SingleJobs/Runtimes/"
    runtimeFile = paste( path, subPath, networkName, "_Cores", cores, "_Runtime.txt", sep = "" )
    if( file.exists( runtimeFile ) )
    {
      write.table( outTable, file = runtimeFile, row.names = FALSE, col.names = FALSE, sep = "\t", dec = ".", quote = FALSE, append = TRUE )
    }
    else
    {
      write.table( outTable, file = runtimeFile, row.names = FALSE, col.names = TRUE, sep = "\t", dec = ".", quote = FALSE, append = FALSE )
    }

    #Show warnings
    if( output )
    {
      warnings()
    }
}
