
#' Data loader module: Create basic folder structure
#'
#' This function creates the standardized folder structure used by regNet for network inference, prediction, and network propagation.
#' @param projectName Name of project folder
#' @param path Directory where the project folder is created
#' @param output Show progress information. Default: TRUE
#' @return String containing the project path.
#' @seealso \code{\link{getProjectPath}}
#' @export
#' @examples 
#' createBasicFolderStructure( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#'
createBasicFolderStructure <- function( projectName, path, output = TRUE )
{
    projectPath = paste( path, projectName, sep = "" )
    
    ##Create project path if it does not exist
    if( !file.exists( projectPath ) )
    {
      dir.create( projectPath )
    }
    else
    {
      if( output )
      {
        print( paste( "Project already exists:", projectPath ) )
        return( projectPath )
      }
    }
    
    ##Create folder structure for network model
    dir.create( paste( projectPath, "/NetworkModel/SingleJobs/Runtimes", sep = "" ), recursive = TRUE )
    dir.create( paste( projectPath, "/NetworkModel/WholeNetwork", sep = "" ) )
    
    ##Create folder for network connectivity
    dir.create( paste( projectPath, "/NetworkConnectivity", sep = "" ) )
    
    ##Create folder for network predictions
    dir.create( paste( projectPath, "/NetworkPredictions", sep = "" ) )
    
    ##Create folder structure for network propagation
    dir.create( paste( projectPath, "/NetworkPropagation/NetworkFlow/BasicNetworkFlowMatrices", sep = "" ), recursive = TRUE )
    dir.create( paste( projectPath, "/NetworkPropagation/NetworkFlow/FinalNetworkFlowMatrices", sep = "" ) )
    dir.create( paste( projectPath, "/NetworkPropagation/ImpactComputations", sep = "" ) )
    
    if( output )
    {
      print( "Create folder structure for project:" )
      print( list.dirs( path = projectPath ) )
    }
    
    return( projectPath )
}


#' Data loader module: Get basic folder structure path
#'
#' This function returns the project path.
#' @param projectName Name of project folder
#' @param path Directory where the project folder is created
#' @return String containing the project path.
#' @export
#' @seealso \code{\link{createBasicFolderStructure}}
#' @examples 
#' getProjectPath( projectName = "MyFirstNetwork", path = "/home/seifert/regNet/AstrocytomaGrades/" )
#'
getProjectPath <- function( projectName, path )
{
    projectPath = paste( path, projectName, sep = "" )
    
    return( projectPath )
}
