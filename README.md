# regNet
regNet is an R package that utilizes gene expression and copy number data to learn regulatory networks for the quantification of potential impacts of individual gene expression alterations on user-defined target genes via network propagation.

# regNet installation
1. Download the folder 'regNet'.
2. Start R and make sure that the libraries 'devtools', 'glmnet', 'Matrix', 'lars', and 'covTest' are installed.
3. Use the following R code to install regNet:

   library( devtools )

   #Set path to the parent directory into which you downloaded the 'regNet' folder
   
   regNetParentDir = "..."
   
   setwd( regNetParentDir )
   
   #
   ### Option 1: Global installation as root into the generally used R package system folder
      
   install( "regNet" )
   
   #
   ### Option 2: Local installation as standard user into a user-specific R package folder     
   #Replace "/home/seifert/LocalRLibs/" in both function calls by your own path
   
   .libPaths( c( .libPaths(), "/home/seifert/LocalRLibs/" ) )
   
   install( pkg = "regNet", args = c( '--library="/home/seifert/LocalRLibs/"' ) )
   
4. regNet should now be installed on your system.

5. Download and unpack the file 'AstrocytomaGrades.zip' from Zenodo at http://doi.org/10.5281/zenodo.580600.
   This file contains the data sets that allow to demonstrate the basic functionality of regNet within 
   a few minutes on a standard computer.

6. Follow the instructions of the R script 'basicCodeUsageExamples.R' to test regNet. See file 'regNet_Vignette.pdf' for
   more details to this case study.
