% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeNetworkFlow.R
\name{computeBasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}
\alias{computeBasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}
\title{Network propagation module: Compute basic initial absolute cohort-specific impact matrix}
\usage{
computeBasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts(data, dataSetName,
  networkName, pValCutoff, localGeneCutoff, path, output = TRUE)
}
\arguments{
\item{data}{Data set containing gene expression and gene copy number profiles for impact computations}

\item{dataSetName}{Name of the data set}

\item{networkName}{Name of the network}

\item{pValCutoff}{Cutoff for significant links. Values from zero (most sparse network) to one (full network) are allowed.}

\item{localGeneCutoff}{Cutoff for removement of regulator links from/to genes in close chromosomal proximity. Set this to zero if no links should be removed.}

\item{path}{Project path}

\item{output}{Show progress information. Default: TRUE}
}
\value{
Saves basic initial impact matrix. See output for details.
}
\description{
This function computes a basic initial cohort-specific impact matrix that is required for the iterative computation of an absolute cohort-specific impact matrix. The basic matrix is computed for a given data set (data) under consideration of a given network (networkName) using network propagation. Only significant regulatory links (pValCutoff) that are not incoming from genes in close chromosomal proximity up- and downstream (localGeneCutoff) of each specific gene are included in the computations. The basic inital impact matrix is saved using pre-defined file naming conventions. This matrix can be used to compute an absolute cohort-specific impact matrix using the function \code{\link{computeNetworkFlow_CohortSpecificAbsoluteImpacts}}. This function allows to separate the computation of the initial from the final impact matrix.
}
\examples{

projectPath = getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
data = loadGeneExpressionAndCopyNumberDataSet( geneExpressionFile = "AS_SignatureTFs_ExpressionLevels.txt", geneCopyNumberFile = "AS_SignatureTFs_CopyNumbers", path = "/home/seifert/regNet/AstrocytomaGrades/Data/" )
computeBasicNetworkFlowMatrix_CohortSpecificAbsoluteImpacts( data = data, dataSetName = "AS_SignatureTFs", networkName = "AS_SignatureTFs", pValCutoff = 0.01, localGeneCutoff = 0, path = projectPath )

}
\seealso{
\code{\link{computeNetworkFlow_CohortSpecificAbsoluteImpacts}}, \code{\link{computeNetworkFlowMatrix_CohortSpecificAbsoluteImpacts}}.
}
