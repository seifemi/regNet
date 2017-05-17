
#' Data loader module: Load data set
#'
#' This function loads a data set (e.g. gene expression, gene copy number) from a given file. The data set must contain the following columns: Gene, Chromosome, Loctcation, Measurements_1, ..., Measurements_N. Gene column: unique gene names. Chromosome column: chromosome entry of genes. Location column: position entry of each gene. Measurement column: gene-specific measurements. Columns in the data set must be separated by tabulators. Genes in the data set must be sorted according to their chromosomal locations.
#' @param file Name of data file
#' @param path Data file path
#' @return List containing measurements (data), genes (genes), chromosomes (chr), and location information (loc).
#' @seealso \code{\link{loadGeneExpressionAndCopyNumberDataSet}}
#' @export
#' @examples
#'
#' data = loadDataSet( loadFile = "AS_SignatureTFs_ExpressionLevels.txt", loadPath = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' 
loadDataSet <- function( file, path )
{
    data = read.delim( file = paste( path, file, sep = "" ), header = TRUE )
    
    M = ncol( data )
    if( M > 4 )
    {
        return( list( data = data[ , 4:M ], genes = as.character( data[ , 1 ] ), chr =  as.character( data[ , 2 ] ), loc = as.numeric( data[ , 3 ] ) ) )
    }
    else
    {
      #Single sample if M == 4
      if( M == 4 )
      {
        return( list( data = matrix( data[ , 4 ], ncol = 1 ), genes = as.character( data[ , 1 ] ), chr =  as.character( data[ , 2 ] ), loc = as.numeric( data[ , 3 ] ) ) )
      }
      else
      {
        return( NULL )
      }
    }
}



#' Data loader module: Load data set that can used by regNet for network inference, prediction, or network propagation
#'
#' This function loads gene expression and corresponding gene copy number data from the given files. The obtained data set can be used by regNet for network learning, prediction, or propagation. See \code{\link{loadDataSet}} for information on how to structure both data files.
#' @param geneExpressionFile Gene expression data file
#' @param geneExpressionFile Gene copy number data file
#' @param path Path to data files
#' @return List containing gene copy number profiles (U), gene expression profiles(Y), genes (genes), chromosomes (chr), and location information (loc).
#' @seealso \code{\link{loadDataSet}} which is called by this function
#' @export
#' @examples 
#'
#' data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#'
loadGeneExpressionAndCopyNumberDataSet <- function( geneExpressionFile, geneCopyNumberFile, path )
{
    ##
    ##Load gene expression and gene copy number data set
    ##
    
    #Expression data
    Y = loadDataSet( file = geneExpressionFile, path = path )
    
    #Copy number data
    U = loadDataSet( file = geneCopyNumberFile, path = path )

    ##  
    ##-Y and U must have the same dimensions and the same order of samples and genes
    ##
    ##Note:
    ## - Please ensure on your own that the genes are sorted in the order of their chromosomal location.
    ## - Please ensure on your own that the order of samples is identical for the gene expression and the gene copy number data set.
    ##
    
    if( ncol( Y$data ) != ncol( U$data ) )
    {
      print( paste( "The gene expression data set '", geneExpressionFile, "' and the gene copy number data set '", geneCopyNumberFile, "' do not contain the same number of samples: ", ncol( Y$data ), " vs. ", ncol( U$data ), sep = "" ) )
      return()
    }

    if( nrow( Y$data ) != nrow( U$data ) )
    {
      print( paste( "The gene expression data set '", geneExpressionFile, "' and the gene copy number data set '", geneCopyNumberFile, "' do not contain the same number of genes: ", nrow( Y$data ), " vs. ", nrow( U$data ), sep = "" ) )
      return()
    }
    
    if( !identical( Y$genes, U$genes ) )
    {
      print( paste( "The gene expression data set '", geneExpressionFile, "' and the gene copy number data set '", geneCopyNumberFile, "' do not contain the genes in the same order.", sep = "" ) )
      return()
    }
    
    if( length( Y$genes ) != length( unique( Y$genes ) ) )
    {
      print( paste( "The gene expression data set '", geneExpressionFile, "' contains the same gene more than once.", sep = "" ) )
      return()
    }
    
    if( length( U$genes ) != length( unique( U$genes ) ) )
    {
      print( paste( "The gene copy number data set '", geneCopyNumberFile, "' contains the same gene more than once.", sep = "" ) )
      return()
    }    
    
    return( list( U = U$data, Y = Y$data, genes = Y$genes, chr = Y$chr, loc = Y$loc ) )
}



#' Data loader module: Test data set compatibility
#'
#' This function evaluates if the two given data sets loaded by \code{\link{loadGeneExpressionAndCopyNumberDataSet}} are compatible. Compatibility of data sets must be ensured to utilize a learned network for predictions or network flow computations.
#' @param testDataSet Test data set
#' @param trainDataSet Data set utilized to learn the network
#' @return Boolean specifying if both data sets are compatible
#' @seealso \code{\link{loadGeneExpressionAndCopyNumberDataSet}}, \code{\link{makeTestDataSetCompatible}}
#' @export
#' @examples
#'
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' testDataSetCompatibility( testDataSet = testData, trainDataSet = trainData )
#'
testDataSetCompatibility <- function( testDataSet, trainDataSet )
{
    ##Test if genes in both data sets are in the same order
    return( identical( testDataSet$genes, trainDataSet$genes ) )
}



#' Data loader module: Transform test data set
#'
#' This function transforms a given test data set to make it compatible to a given training data set that was used to learn a network. Additional genes are removed, missing genes are artificially added with measurement values of zero, and the same order of genes is ensured.
#' @param testDataSet Test data set
#' @param trainDataSet Data set used for network inference
#' @param output Show progress information. Default: TRUE
#' @return Test data set compatible with training data set
#' @seealso \code{\link{loadGeneExpressionAndCopyNumberDataSet}}, \code{\link{testDataSetCompatibility}}
#' @export
#' @examples
#'
#' testData  = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt", geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' trainData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
#' gbmTFData = makeTestDataSetCompatible( testDataSet = testData, trainDataSet = trainData )
#'
makeTestDataSetCompatible <- function( testDataSet, trainDataSet, output = TRUE )
{
    #Both data sets arre compatible
    if( testDataSetCompatibility( testDataSet = testDataSet, trainDataSet = trainDataSet ) )
    {
      return( testDataSet )
    }
    
    
    ##
    ##Transform test data sets
    ##
    genesTestDataSet  = testDataSet$genes
    genesTrainDataSet = trainDataSet$genes
    
    ##Remove additional genes
    additionalGenesInTestDataSet = setdiff( genesTestDataSet, genesTrainDataSet )
    if( length( additionalGenesInTestDataSet ) > 0 )
    {
      #Determine positions of additional genes
      pos = match( additionalGenesInTestDataSet, genesTestDataSet )

      #Remove additional genes
      testDataSet$U     = testDataSet$U[ -pos, ]
      testDataSet$Y     = testDataSet$Y[ -pos, ]
      testDataSet$genes = testDataSet$genes[ -pos ]
      testDataSet$chr   = testDataSet$chr[ -pos ]
      testDataSet$loc   = testDataSet$loc[ -pos ]

      if( output )
      {
        print( "Genes removed from test data set:" )
        print( additionalGenesInTestDataSet )
      }
    }
    
    ##Add missing genes
    missingGenesInTestDataSet = setdiff( genesTrainDataSet, genesTestDataSet )
    if( length( missingGenesInTestDataSet ) > 0 )
    {
      #Determine positions of missing genes
      pos = match( missingGenesInTestDataSet, genesTrainDataSet )

      #Add missing genes to test data set
      dummyMat = matrix( 0, ncol = ncol( testDataSet$U ), nrow = length( pos ) ) #U and Y have same number of columns and rows
      colnames( dummyMat ) <- colnames( testDataSet$U )
      testDataSet$U     = rbind( testDataSet$U, dummyMat ) ##Add missing genes with artificial measurement values of zero
      colnames( dummyMat ) <- colnames( testDataSet$Y )
      testDataSet$Y     = rbind( testDataSet$Y, dummyMat )
      testDataSet$genes = c( testDataSet$genes, trainDataSet$genes[ pos ] )
      testDataSet$chr   = c( testDataSet$chr, trainDataSet$chr[ pos ] )
      testDataSet$loc   = c( testDataSet$loc, trainDataSet$loc[ pos ] )

      if( output )
      {
        print( "Genes added to test data set:" )
        print( missingGenesInTestDataSet )
      }
    }
    
    ##Ensure same order of genes in both data sets
    if( output )
    {
      print( "Reorder genes in test data set." )
    }
    
    genesTestDataSet  = testDataSet$genes
    genesTrainDataSet = trainDataSet$genes
    
    pos = match( genesTrainDataSet, genesTestDataSet )
    if( !is.vector( testDataSet$U ) )
    {
      #More than one sample: U and Y are matrices
      testDataSet$U     = testDataSet$U[ pos, ]
      testDataSet$Y     = testDataSet$Y[ pos, ]
    }
    else
    {
      #Only one sample: U and Y are vectors that must be transformed back to matrices
      L = length( testDataSet$U )
      testDataSet$U   = as.matrix( testDataSet$U[ pos ], nrow = L, ncol = 1 )
      testDataSet$Y   = as.matrix(  testDataSet$Y[ pos ], nrow = L, ncol = 1 )
    }
    testDataSet$genes = testDataSet$genes[ pos ]
    testDataSet$chr   = testDataSet$chr[ pos ]
    testDataSet$loc   = testDataSet$loc[ pos ]

    return( testDataSet )
}
