## First part, install required libraries
if( !require( randomForest ) ) install.packages("randomForest");require(randomForest)
if( !require( randomForestSRC ) ) install.packages("randomForestSRC");require(randomForestSRC)
if( !require( parallel ) ) install.packages("parallel");require(parallel)
if( !require( survival ) ) install.packages("survival");require(survival)
if( !require( gatbxr ) ){
  if( !require( devtools ) ) install.packages("devtools")
  devtools::install_github('drizztxx/gatbxr')
  require(gatbxr)
}
## Second part, initialized siga environment
if('.SIGA' %in% search()) detach('.SIGA')
if(exists('.SIGA')) rm(.SIGA)
if(exists('.tempVars')) rm(.tempVars)
.SIGA = new.env()
.tempVars = new.env()

assign('pathIntrospection', function(){
  frameFiles = Filter(Negate(is.null),
                      lapply(sys.frames(),function(x)x$ofile))
  return(dirname(frameFiles[[length(frameFiles)]]))},
  envir = .tempVars
)
.tempVars$loadpath = .tempVars$pathIntrospection()

.tempVars$rfiles = c("SIGA.tte.R","VariableScreen.R","parSIGA.R","permadjp.R",
                     "predict.SIGA.R","SIGA.binary.R","SIGA_CV.R","Best.R",
                     "decodsubset.R","objectiveFunctions.R","boot.SIGA.R","freqvars.R",
                     "SIGA.R")
.tempVars$path2rfiles = sapply(.tempVars$rfiles,function(x){
  ifelse(.tempVars$loadpath!='',file.path(.tempVars$loadpath,x),x)})

sapply(.tempVars$path2rfiles,function(x)sys.source(x,.SIGA))
attach(.SIGA,warn.conflicts=FALSE)
rm(.tempVars)

