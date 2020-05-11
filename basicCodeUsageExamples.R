#######################################################################################################################################################################
##
##Basic regNet code usage examples.
##
#######################################################################################################################################################################
##
##Contact: michael.seifert@tu-dresden.de
##


##
##Set library path to installed regNet package and load package
##

##You have to set your specific path! (See https://github.com/seifemi/regNet for details how to install regNet.)
libLoc = "/home/seifert/Documents/Latex/PaperProjects/regNet_Rpackage/regNet_R_package/regNetTestInstallation_LocalLibs/"
library( regNet, lib.loc = libLoc )


##
##Set your path to the project
##

##You have to set your specific path! (Data of file 'AstrocytomaGrades.zip' from Zenodo at http://doi.org/10.5281/zenodo.580600.)
myPath = "/home/seifert/AstrocytomaGrades/"


##Fixed random seed
#set.seed( 123 )


##
##Print standard output
##
output = TRUE

#######################################################################################################################################################################
##
## Create basic folder structure
##
#######################################################################################################################################################################
projectName = "MyFirstNetwork"
path = myPath

projectPath = createBasicFolderStructure( projectName = projectName, path = path, output = output )


 
####################################################################################################################################################################### 
##
## Load data set
##
#######################################################################################################################################################################
geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt"
geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers.txt"
loadPath = paste0( myPath, "Data/" )
  
data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = geneExpressionFile, geneCopyNumberFile = geneCopyNumberFile, path = loadPath )



#######################################################################################################################################################################
##
## Network inference
##
#######################################################################################################################################################################
networkName = "AS_SignatureTFs"
totalNumberOfJobs = 4

for( i in 1:totalNumberOfJobs )
{
    learnNetwork_ParallelComputation( data = data, networkName = networkName, cores = totalNumberOfJobs, job = i, path = projectPath, nfolds = 10, cvReplicates = 10, output = output )
}


##Combine single jobs
combineSingleJobs( networkName = networkName, cores = totalNumberOfJobs, path = projectPath, output = output )

##
##Get network connectivity table
##
getNetworkConnections( networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath, output = output )

##
##Create random network instances
##
for( i in 1:10 )
{
    determineRandomNetworkWithFilteringForSignificantPredictors( networkName = networkName, pValCutoff = 0.01, randomNetworkInstance = i, path = projectPath, output = output )
}

#######################################################################################################################################################################



#######################################################################################################################################################################
##
## Network-based predictions of gene expression levels
##
#######################################################################################################################################################################

##
##TCGA GBM test data
##
geneExpressionFile = "TCGA_GBM_ExpressionLevels.txt"
geneCopyNumberFile = "TCGA_GBM_CopyNumbers.txt"
loadPath = paste0( myPath, "Data/" )

gbmData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = geneExpressionFile, geneCopyNumberFile = geneCopyNumberFile, path = loadPath )
print( testDataSetCompatibility( testDataSet = gbmData, trainDataSet = data ) )

gbmTFData = makeTestDataSetCompatible( testDataSet = gbmData, trainDataSet = data )
print( testDataSetCompatibility( testDataSet = gbmTFData, trainDataSet = data ) )

networkName = "AS_SignatureTFs"
dataSetName = "TCGA_GBM_SignatureTFs"
w <- getOption( "warn" )
options( warn = -1 )
predictGeneExpression( data = gbmTFData, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath, output = output )
options( warn = w )

##Random network instances on GBM data
for( i in 1:10 )
{
    randomNetworkName = paste0( "RandomNetwork_", i, "_PValueCutoff_0.01_BasedOn_", networkName )
    dataSetName = "TCGA_GBM_SignatureTFs"
    predictGeneExpression( data = gbmTFData, dataSetName = dataSetName, networkName = randomNetworkName, pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath, output = output )
}


##
##TCGA LGG test data
##
geneExpressionFile = "TCGA_LGG_ExpressionLevels.txt"
geneCopyNumberFile = "TCGA_LGG_CopyNumbers.txt"
loadPath = paste0( myPath, "Data/" )

lggData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = geneExpressionFile, geneCopyNumberFile = geneCopyNumberFile, path = loadPath )
print( testDataSetCompatibility( testDataSet = lggData, trainDataSet = data ) )

lggTFData = makeTestDataSetCompatible( testDataSet = lggData, trainDataSet = data )
print( testDataSetCompatibility( testDataSet = lggTFData, trainDataSet = data ) )

networkName = "AS_SignatureTFs"
dataSetName = "TCGA_LGG_SignatureTFs"
w <- getOption( "warn" )
options( warn = -1 )
predictGeneExpression( data = lggTFData, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath, output = output )
options( warn = w )

##Random network instances on LGG data
for( i in 1:10 )
{
    randomNetworkName = paste0( "RandomNetwork_", i, "_PValueCutoff_0.01_BasedOn_", networkName )
    dataSetName = "TCGA_LGG_SignatureTFs"
    predictGeneExpression( data = lggTFData, dataSetName = dataSetName, networkName = randomNetworkName, pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath, output = output )
}


##
##PA test data
##
geneExpressionFile = "PA_GSE5675_ExpressionLevels.txt"
geneCopyNumberFile = "PA_GSE5675_CopyNumbers.txt"
loadPath = paste0( myPath, "Data/" )

paData = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = geneExpressionFile, geneCopyNumberFile = geneCopyNumberFile, path = loadPath )
print( testDataSetCompatibility( testDataSet = paData, trainDataSet = data ) )

paTFData = makeTestDataSetCompatible( testDataSet = paData, trainDataSet = data )
print( testDataSetCompatibility( testDataSet = paTFData, trainDataSet = data ) )

networkName = "AS_SignatureTFs"
dataSetName = "PA_GSE5675_SignatureTFs"
w <- getOption( "warn" )
options( warn = -1 )
predictGeneExpression( data = paTFData, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath, output = output )
options( warn = w )

##Random network instances on PA data
for( i in 1:10 )
{
    randomNetworkName = paste0( "RandomNetwork_", i, "_PValueCutoff_0.01_BasedOn_", networkName )
    dataSetName = "PA_GSE5675_SignatureTFs"
    predictGeneExpression( data = paTFData, dataSetName = dataSetName, networkName = randomNetworkName, pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath, output = output )
}


##
##Prediction summary plot
##
x11()
testDataSets = c( "TCGA_GBM_SignatureTFs", "TCGA_LGG_SignatureTFs", "PA_GSE5675_SignatureTFs" )
titles = c( "GBM: Prediction quality", "LGG: Prediction quality", "PA: Prediction quality" )
layout( mat = matrix( 1:4, nrow = 2, ncol = 2 ) )

for( i in 1:3 )
{
    orig_pred = as.numeric( read.delim( file = paste( projectPath, "/NetworkPredictions/", testDataSets[ i ], "_PredictionOfGeneExpressionBasedOn_AS_SignatureTFs_PValueCutoff_0.01_LocalGeneCutoff_0.txt", sep = "" ), header = TRUE )[ , 2 ] )
    
    avgRandPredictions = c()
    for( r in 1:10 )
    {
      randPred = as.numeric( read.delim( file = paste( projectPath, "/NetworkPredictions/", testDataSets[ i ], "_PredictionOfGeneExpressionBasedOn_RandomNetwork_", i, "_PValueCutoff_0.01_BasedOn_AS_SignatureTFs_PValueCutoff_0.01_LocalGeneCutoff_0.txt", sep = "" ), header = TRUE )[ , 2 ] )
      if( r == 1 )
      {
        avgRandPredictions = randPred
      }
      else
      {
        avgRandPredictions = avgRandPredictions + randPred
      }
    }
    rand_pred = avgRandPredictions / 10
    
    hist( rand_pred, breaks = seq( -1, 1, length.out = 40 ), xlim = c( -1, 1 ), xlab = "Correlation", col = "darkgrey", main = titles[ i ] )
    hist( orig_pred, breaks = seq( -1, 1, length.out = 40 ), xlim = c( -1, 1 ), add = TRUE, col = rgb( red = 0.85, green = 0, blue = 0, alpha = 0.7 ) )
    legend( x = -1.2, y = 10, legend = c( "TF-Network", "Random" ), text.col = c( rgb( red = 0.85, green = 0, blue = 0, alpha = 0.7 ), "darkgrey" ), lty = 0, col = c( rgb( red = 0.85, green = 0, blue = 0, alpha = 0.7 ), "darkgrey" ), bty = "n" )
}



#######################################################################################################################################################################
##
## Network propagation: Compute cohort-specific impacts
##
#######################################################################################################################################################################

##
##
##Which hub TF has the greatest impact on all other reachable TFs
##
##

##
##Get average impacts: cohort-specific absolute impacts
##
dataSetName = "AS_SignatureTFs"
networkName = "AS_SignatureTFs"

##Compute cohort-specific impact matrix for all gene pairs for the given data set (data: paired gene expression and copy number profiles) with respect to the learned network (networkName)
computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts( data = data, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, output = output )

##Get average impact of each source gene (sourceGenes) on all target genes (targetGenes) from the previously computed cohort-specific impact matrix
sourceGenes = c( "RBBP4", "NFIA", "MEOX2", "PAX6", "ZNF337", "THRB", "ZCCHC24", "CCNL2", "TBR1", "ZNF300", "APBA1", "GPR123" )
targetGenes = data$genes
outputFile = "HubTFs_AverageAbsoluteImpactsOnOtherGenes.txt"
res_orig = getAverageImpacts_CohortSpecificAbsoluteImpacts( sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

##Get impacts of each source gene (sourceGenes) on each target gene (targetGenes) from the previously computed cohort-specific impact matrix
outputFile = "HubTFs_AbsoluteImpactsOnOtherGenes.txt"
res = getImpacts_CohortSpecificAbsoluteImpacts( sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

##Get corresponding average impacts under ten random networks
avgRandImpact = c()
for( i in 1:10 )
{
    randomNetworkName = paste0( "RandomNetwork_", i, "_PValueCutoff_0.01_BasedOn_", networkName )
    dataSetName = "AS_SignatureTFs"
    
    ##Compute cohort-specific impact matrix for all gene pairs for the given data set with respect to the random network (randomNetworkName)
    computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts( data = data, dataSetName = dataSetName, networkName = randomNetworkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, output = output )
    
    ##Get average impact of each source gene on all target genes from the cohort-specific impact matrix previously computed for the random network
    outputFile = paste0( "HubTFs_AverageAbsoluteImpactsOnOtherGenesRandomNetwork_", i, ".txt" )
    dummy = getAverageImpacts_CohortSpecificAbsoluteImpacts( sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = randomNetworkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )
    
    if( i == 1 )
    {
      avgRandImpact = dummy[ , 2 ]
    }
    else
    {
      avgRandImpact = avgRandImpact + dummy[ , 2 ]
    }
}
avgRandImpact = avgRandImpact / 10


##
##Impact summary plot for average impacts
##
L = length( res_orig[ , 2 ] )
plot( x = 1:L, y = res_orig[ , 2 ], type = "h", col = rgb( red = 0.85, green = 0, blue = 0, alpha = 0.7 ), xlab = "", ylab = "Average impact", main = "Hub-regulator impacts", axes = FALSE )        
points( x = 1:L + 0.2, y = avgRandImpact, type = "h", col = "darkgrey" )
legend( x = 0, y = 0.08, legend = c( "TF-Network", "Random" ), text.col = c( rgb( red = 0.85, green = 0, blue = 0, alpha = 0.7 ), "darkgrey" ), lty = 0, col = c( rgb( red = 0.85, green = 0, blue = 0, alpha = 0.7 ), "darkgrey" ), bty = "n" )
axis( 1, at = 1:L, labels = sourceGenes, las = 2, cex.axis = 0.8 )
axis( 2 )



#######################################################################################################################################################################
##
## Network propagation: Compute patient-specific absolute impacts
##
#######################################################################################################################################################################

##
##
##Which hub TF has the greatest impact on all other reachable TFs for selected patients
##
##

##
##Get average impacts: patient-specific absolute impacts
##
networkName = "AS_SignatureTFs"

##
##Compute patient-specific impact matrix for all gene pairs for a specific patient in the given data set with respect to the previously learned network
##

##Impact matrix for patient 6 of the TCGA GBM TF data set
dataSetName = "TCGA_GBM_SignatureTFs"
w <- getOption( "warn" )
options( warn = -1 )
computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts( patient = 6, data = gbmTFData, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, output = output )
options( warn = w )

##Get average absolute impact of each source gene (sourceGenes) on all target genes (targetGenes) from the previously computed patient-specific impact matrix for patient 6 of the TCGA GBM TF data set
sourceGenes = c( "RBBP4", "NFIA", "MEOX2", "PAX6", "ZNF337", "THRB", "ZCCHC24", "CCNL2", "TBR1", "ZNF300", "APBA1", "GPR123" )
targetGenes = data$genes
outputFile = "HubTFs_AverageAbsoluteImpactsOnOtherGenes_Patient_6_TCGA_GBM_DataSet.txt"
resPatient6_GBM_abs = getAverageImpacts_PatientSpecificAbsoluteImpacts( patient = 6, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

##Get absolute impacts of each source gene (sourceGenes) on each target gene (targetGenes) from the previously computed patient-specific impact matrix for patient 6 of the TCGA GBM TF data set
outputFile = "HubTFs_AbsoluteImpactsOnOtherGenes_Patient_6_TCGA_GBM_DataSet.txt"
res = getImpacts_PatientSpecificAbsoluteImpacts( patient = 6, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )


##Impact matrix for patient 18 of the PA TF data set
dataSetName = "PA_GSE5675_SignatureTFs"
computeNetworkFlowMatrix_PatientSpecificAbsoluteImpacts( patient = 18, data = paTFData, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, output = output )

##Get average absolute impact of each source gene (sourceGenes) on all target genes (targetGenes) from the previously computed patient-specific impact matrix for patient 18 of the PA TF data set
sourceGenes = c( "RBBP4", "NFIA", "MEOX2", "PAX6", "ZNF337", "THRB", "ZCCHC24", "CCNL2", "TBR1", "ZNF300", "APBA1", "GPR123" )
targetGenes = data$genes
outputFile = "HubTFs_AverageAbsoluteImpactsOnOtherGenes_Patient_18_PA_DataSet.txt"
resPatient18_PA_abs = getAverageImpacts_PatientSpecificAbsoluteImpacts( patient = 18, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

##Get absolute impacts of each source gene (sourceGenes) on each target gene (targetGenes) from the previously computed patient-specific impact matrix for patient 18 of the PA TF data set
outputFile = "HubTFs_AbsoluteImpactsOnOtherGenes_Patient_18_PA_DataSet.txt"
res = getImpacts_PatientSpecificAbsoluteImpacts( patient = 18, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )


##
##Impact summary plot for average impacts
##
x11( width = 11, height = 5 )
layout( mat = matrix( 1:2, nrow = 1, ncol = 2 ) )
L = length( resPatient6_GBM_abs[ , 2 ] )
plot( x = 1:L, y = resPatient6_GBM_abs[ , 2 ], type = "h", col = "red", xlab = "", ylab = "Average absolute impact", main = "Patient-specific absolute hub impacts", axes = FALSE )        
points( x = 1:L + 0.2, y = resPatient18_PA_abs[ , 2 ], type = "h", col = "green" )
legend( x = 0, y = 0.035, legend = c( "GBM: Patient 6", "PA: Patient 18" ), text.col = c( "red", "green" ), lty = 0, col = c( "red", "green" ), bty = "n" )
axis( 1, at = 1:L, labels = sourceGenes, las = 2, cex.axis = 0.8 )
axis( 2 )



#######################################################################################################################################################################
##
## Network propagation: Compute patient-specific relative impacts
##
#######################################################################################################################################################################

##
##
##Which hub TFs act on average as activators (positive impact) or inhibitors (negative impact) of all other reachable TFs for selected patients
##
##

##
##Get average impacts: patient-specific relative impacts
##
networkName = "AS_SignatureTFs"

##
##Compute patient-specific impact matrix for all gene pairs for a specific patient in the given data set with respect to the previously learned network
##

##Impact matrix for patient 6 of the TCGA GBM TF data set
dataSetName = "TCGA_GBM_SignatureTFs"
w <- getOption( "warn" )
options( warn = -1 )
computeNetworkFlowMatrix_PatientSpecificRelativeImpacts( patient = 6, data = gbmTFData, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, output = output )
options( warn = w )

##Get average relative impact of each source gene (sourceGenes) on all target genes (targetGenes) from the previously computed patient-specific impact matrix for patient 6 of the TCGA GBM TF data set
sourceGenes = c( "RBBP4", "NFIA", "MEOX2", "PAX6", "ZNF337", "THRB", "ZCCHC24", "CCNL2", "TBR1", "ZNF300", "APBA1", "GPR123" )
targetGenes = data$genes
outputFile = "HubTFs_AverageRelativeImpactsOnOtherGenes_Patient_6_TCGA_GBM_DataSet.txt"
resPatient6_GBM_rel = getAverageImpacts_PatientSpecificRelativeImpacts( patient = 6, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

##Get relative impacts of each source gene (sourceGenes) on each target gene (targetGenes) from the previously computed patient-specific impact matrix for patient 6 of the TCGA GBM TF data set
outputFile = "HubTFs_RelativeImpactsOnOtherGenes_Patient_6_TCGA_GBM_DataSet.txt"
res = getImpacts_PatientSpecificRelativeImpacts( patient = 6, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )


##Impact matrix for patient 18 of the PA TF data set
dataSetName = "PA_GSE5675_SignatureTFs"
computeNetworkFlowMatrix_PatientSpecificRelativeImpacts( patient = 18, data = paTFData, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, output = output )

##Get average relative impact of each source gene (sourceGenes) on all target genes (targetGenes) from the previously computed patient-specific impact matrix for patient 18 of the PA TF data set
sourceGenes = c( "RBBP4", "NFIA", "MEOX2", "PAX6", "ZNF337", "THRB", "ZCCHC24", "CCNL2", "TBR1", "ZNF300", "APBA1", "GPR123" )
targetGenes = data$genes
outputFile = "HubTFs_AverageRelativeImpactsOnOtherGenes_Patient_18_PA_DataSet.txt"
resPatient18_PA_rel = getAverageImpacts_PatientSpecificRelativeImpacts( patient = 18, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

##Get relative impacts of each source gene (sourceGenes) on each target gene (targetGenes) from the previously computed patient-specific impact matrix for patient 18 of the PA TF data set
outputFile = "HubTFs_RelativeImpactsOnOtherGenes_Patient_18_PA_DataSet.txt"
res = getImpacts_PatientSpecificRelativeImpacts( patient = 18, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )


##
##Impact summary plot for average impacts
##
L = length( resPatient6_GBM_rel[ , 2 ] )
plot( x = 1:L, y = resPatient6_GBM_rel[ , 2 ], type = "h", col = "red", xlab = "", ylab = "Average relative impact", main = "Patient-specific relative hub impacts", axes = FALSE )        
points( x = 1:L + 0.2, y = resPatient18_PA_rel[ , 2 ], type = "h", col = "green" )
lines( x = c( 0, L ), y = c( 0, 0 ), lty = 3, lwd = 0.3 )
legend( x = 0, y = 0.02, legend = c( "GBM: Patient 6", "PA: Patient 18" ), text.col = c( "red", "green" ), lty = 0, col = c( "red", "green" ), bty = "n" )
axis( 1, at = 1:L, labels = sourceGenes, las = 2, cex.axis = 0.8 )
axis( 2 )




#######################################################################################################################################################################
##
## Network propagation: Compute impacts for a new cohort or a new patient utilizing a pre-computed correlation statistics
##
#######################################################################################################################################################################

##
##
##Let us consider the following situation:
##- We have already learned a regulatory network from a large cohort of patients.
##- We evaluated the predictive power of the learned network by predicting the expression levels of genes for another large test cohort of patients with similar
##  gene expression and gene copy number profiles (e.g. patients that were diagnosed to have a specific type of cancer). Thus, we already have the correlation
##  statistics for the test cohort containing the gene-specific correlations between predicted and originally measured expression levels.
##- Now, a new little cohort or a single patient has to be analyzed, but the number of samples in this new data set is not large enough to obtain a robust correlation
##  statistics for these new data.
##- regNet allows to analyze this new data sets by utilizing a learned network and a pre-computed correlation statistics.
##
##The following code examples demonstrate how this can be used to obtain cohort-specific and patient-specific impact scores for a new cohort or patient.
##

######################
##
##
##Cohort-specific impacts for new data set
##
##
networkName = "AS_SignatureTFs"
corStatDataSetName = "TCGA_LGG_SignatureTFs"

##
##Compute patient-specific absolute impact matrix for all gene pairs for a given data set (gbmTFData) with respect to the previously learned network (networkName) and the previously computed correlation statistics (corStatDataSetName)
##

##Impact matrix for TCGA GBM TF data set
dataSetName = "TCGA_GBM_SignatureTFs"
w <- getOption( "warn" )
options( warn = -1 )
computeNetworkFlowMatrixBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( data = gbmTFData, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, output = output )
options( warn = w )

##Get average absolute impact of each source gene (sourceGenes) on all target genes (targetGenes) from the previously computed impact matrix
sourceGenes = c( "RBBP4", "NFIA", "MEOX2", "PAX6", "ZNF337", "THRB", "ZCCHC24", "CCNL2", "TBR1", "ZNF300", "APBA1", "GPR123" )
targetGenes = data$genes
outputFile = "HubTFs_AverageAbsoluteImpactsOnOtherGenes_Utilizing_CorStat_TCGA_LGG_SignatureTFs.txt"
resPatient6_GBM_rel = getAverageImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

##Get absolute impacts of each source gene (sourceGenes) on each target gene (targetGenes) from the previously computed impact matrix
outputFile = "HubTFs_AbsoluteImpactsOnOtherGenes_Utilizing_CorStat_TCGA_LGG_SignatureTFs.txt"

res = getImpactsBasedOnGivenCorStat_CohortSpecificAbsoluteImpacts( sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, corStatDataSetName = corStatDataSetName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

######################


######################
##
##
##Patient-specific absolute impacts for a new patient
##
##
networkName = "AS_SignatureTFs"
corStatDataSetName = "TCGA_LGG_SignatureTFs"

##
##Compute patient-specific absolute impact matrix for all gene pairs for a specific patient in the given data set (gbmTFData) with respect to the previously learned network (networkName) and the previously computed correlation statistics (corStatDataSetName)
##

##Impact matrix for patient 6 of the TCGA GBM TF data set
dataSetName = "TCGA_GBM_SignatureTFs"
w <- getOption( "warn" )
options( warn = -1 )
computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = 6, data = gbmTFData, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, output = output )
options( warn = w )

##Get average absolute impact of each source gene (sourceGenes) on all target genes (targetGenes) from the previously computed patient-specific impact matrix for patient 6 of the TCGA GBM TF data set
sourceGenes = c( "RBBP4", "NFIA", "MEOX2", "PAX6", "ZNF337", "THRB", "ZCCHC24", "CCNL2", "TBR1", "ZNF300", "APBA1", "GPR123" )
targetGenes = data$genes
outputFile = "HubTFs_AverageAbsoluteImpactsOnOtherGenes_Patient_6_TCGA_GBM_DataSet_Utilizing_CorStat_TCGA_LGG_SignatureTFs.txt"
resPatient6_GBM_rel = getAverageImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = 6, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

##Get absolute impacts of each source gene (sourceGenes) on each target gene (targetGenes) from the previously computed patient-specific impact matrix for patient 6 of the TCGA GBM TF data set
outputFile = "HubTFs_AbsoluteImpactsOnOtherGenes_Patient_6_TCGA_GBM_DataSet_Utilizing_CorStat_TCGA_LGG_SignatureTFs.txt"

res = getImpactsBasedOnGivenCorStat_PatientSpecificAbsoluteImpacts( patient = 6, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, corStatDataSetName = corStatDataSetName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

######################


######################
##
##
##Patient-specific relative impacts for a new patient
##
##
networkName = "AS_SignatureTFs"
corStatDataSetName = "TCGA_LGG_SignatureTFs"

##
##Compute patient-specific relative impact matrix for all gene pairs for a specific patient the given data set (gbmTFData) with respect to the previously learned network (networkName) and the previously computed correlation statistics (corStatDataSetName)
##

##Impact matrix for patient 6 of the TCGA GBM TF data set
dataSetName = "TCGA_GBM_SignatureTFs"
w <- getOption( "warn" )
options( warn = -1 )
computeNetworkFlowMatrixBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = 6, data = gbmTFData, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, output = output )
options( warn = w )

##Get average relative impact of each source gene (sourceGenes) on all target genes (targetGenes) from the previously computed patient-specific impact matrix for patient 6 of the TCGA GBM TF data set
sourceGenes = c( "RBBP4", "NFIA", "MEOX2", "PAX6", "ZNF337", "THRB", "ZCCHC24", "CCNL2", "TBR1", "ZNF300", "APBA1", "GPR123" )
targetGenes = data$genes
outputFile = "HubTFs_AverageRelativeImpactsOnOtherGenes_Patient_6_TCGA_GBM_DataSet_Utilizing_CorStat_TCGA_LGG_SignatureTFs.txt"
resPatient6_GBM_rel = getAverageImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = 6, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, corStatDataSetName = corStatDataSetName, networkName = networkName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

##Get relative impacts of each source gene (sourceGenes) on each target gene (targetGenes) from the previously computed patient-specific impact matrix for patient 6 of the TCGA GBM TF data set
outputFile = "HubTFs_RelativeImpactsOnOtherGenes_Patient_6_TCGA_GBM_DataSet_Utilizing_CorStat_TCGA_LGG_SignatureTFs.txt"

res = getImpactsBasedOnGivenCorStat_PatientSpecificRelativeImpacts( patient = 6, sourceGenes = sourceGenes, targetGenes = targetGenes, dataSetName = dataSetName, networkName = networkName, corStatDataSetName = corStatDataSetName, pValCutoff = 0.01, localGeneCutoff = 0, colSumsThreshold = 1e-3, path = projectPath, outputFile = outputFile, output = output )

######################

