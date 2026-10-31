SIGA.binary <- function(data,
                        trt="treat=='E'",
                        ctrl="treat=='C'",
                        resp="res==1",
                        varScreen=c(1,5),
                        varused=NULL,
                        popSize=100,
                        subpop=1,
                        objFunction=c("objF.glm2","objF.glm","objF.difference","objF.difference2"),
                        MAXGEN = 400,        ## max Number of generations
                        GGAP = 0.9,          ## Generation gap, how many new individuals are created
                        SEL_F = 'sus',       ## Name of selection function
                        XOV_F = 'reclin',     ## Name of recombination function for individuals
                        MUT_F = 'mutbga',       ## Name of mutation function for individuals
                        cutoff = 10,
                        verbose = TRUE,
                        ...
){
  library(randomForest)
  trtvar <- strsplit(trt,'=')[[1]][1]
  resvar <- strsplit(resp,'=')[[1]][1]
  if (!is.null(varused)){
    col.vim <- varused
    top <- NROW(col.vim)
  } else if (varScreen[1]) { ## Screen the top important variables
    top <- varScreen[2]
    screen.res <- VariableScreen(data=data,trt=trt,ctrl=ctrl,resp=resp,varScreen = varScreen,resptyp = "binary")
    col.vim <- screen.res$col.vim[1:top]
  } else {
    top <- NCOL(data) - 2
    col.vim <- colnames(data)[!colnames(data) %in% c(trtvar,resvar)]
  }
  
  
  
  
  ## Create the initial population used crtrp method.
  ## Each chrom of the initial population stands for a possible solution(subgroup), and will have
  ## #top# genes. Each gene will represent the status of each selected top variable.
  ## The genetic representation rule is as below:
  ## 1. Value of gene(B) is ranged from 0-4 (uniform distribution).
  ## 2. For continuous predictor, if B <= 2, it will not be selected (50% possibility);
  ##    if 2 < B <=3, then the subgroup will be defined as <= (B - 2) quantile;
  ##    if 3 < B <=4, then the subgroup will be defined as >= (B - 3) quantile.
  ## 3. For binary predictor, B <= 2 indicates the variable will not be selected;
  ##    and 2 < B <= 3 indicates the first category(level of factor) will be used to define the subgroup 
  ##    whilst 3 < B <= 4 indicates the second one will be used.
  ## 4. For multi-categorical predictor (i=1 to C), it will not be selected if B <= 2;
  ##    else there will be (2^C - 2) subgroups to be choosed and each subgroup is one possible
  ##    combination of the category for the predictor. The value of B (2-4) will be equally divided to
  ##    (2^C - 2) blocks and each block stands for one category combination.
  FieldDR <- matrix(c(0,4),2,top)
  Chrom <- crtrp(popSize*subpop,FieldDR)
  Chrom.pool <- list()
  Chrom.pool[[1]] <- Chrom
  
  ## Decoding the generated Chrom according to the genetic representation rule.
  decod.subset <- decodsubset(data=data[,col.vim],Chrom = Chrom)
  subset.pool <- list()
  subset.pool[[1]] <- decod.subset$subgroup.desc
  
  ## Caculate the treatment effect (that is the pvalue caculated through the glm model) for every subgroups
   ObjV <- do.call(objFunction,list(subgroups=decod.subset$subgroup.desc, 
                                   data=data[,c(col.vim,resvar,trtvar)],
                                   cutoff,trt,ctrl,resp,resvar,trtvar,... ))
  ##fun_wrapper <- function(subgroups,e=parent.frame()){
  ##  lapply(objFunction,do.call,args=list(subgroups=subgroups, 
  ##                                       data=data[,c(col.vim,resvar,trtvar)],trt = trt,
  ##                                       ctrl = ctrl,resp = resp, ... ),envir=e)
  ##}
  ##ObjV <- fun_wrapper(decod.subset$subgroup.desc)
  
  objv.pool <- list()
  objv.pool[[1]] <- ObjV$objv
  objv <- ObjV$objv
  
  ## reset count variables
  gen = 0;
  Best = rep(NA,MAXGEN)
  
  ## Iterate population
  while (gen < MAXGEN){
    ## Calculate objective function for population
    Best[gen+1] = min(objv);
    if (verbose) cat(sprintf("Best objective value of the %s generation is %s \n",gen+1,Best[gen+1]))
    
    ## Fitness assignement to whole population
    FitnV = ranking(objv,SUBPOP = subpop);
    
    ## Select individuals from population
    SelCh = gatbxr::select(SEL_F, Chrom, FitnV, GGAP, subpop);
    
    ## Recombine selected individuals (crossover)
    SelCh=recombin("reclin",SelCh,SUBPOP = subpop)
    
    ## Mutate offspring
    SelCh=gatbxr::mutate(MUT_F, SelCh, FieldDR = FieldDR, SUBPOP = subpop);
    
    ##Evaluate offspring, call objective function
    ObjVSel <- do.call(objFunction,list(subgroups=decodsubset(data=data[,col.vim],Chrom = SelCh)$subgroup.desc, 
                                        data=data[,c(col.vim,trtvar,resvar)],
                                        cutoff,trt,ctrl,resp,resvar,trtvar,... ))
    ##ObjVSel <- fun_wrapper(decodsubset(data=data[,col.vim],Chrom = SelCh)$subgroup.desc)
    
    
    ## Insert offspring in population replacing parents
    rs = reins(Chrom, SelCh,subpop,c(1,1),objv,ObjVSel$objv);
    Chrom = rs$Chrom;
    objv  = rs$ObjVCh
    
    gen = gen + 1;
    Chrom.pool[[gen]] <- Chrom
    objv.pool[[gen]] <- objv
    subset.pool[[gen]] <- decodsubset(data=data[,col.vim],Chrom = Chrom)$subgroup.desc
  }
  ## End of script
  dsn <- deparse(substitute(data))
  if (!exists("screen.res")) screen.res <- NULL
  final.res <- list(Chrom=Chrom.pool,Objv=objv.pool,Subset=subset.pool,Best=Best,varused=col.vim,dsn=dsn,resptyp="binary",
                    screen.res=screen.res,call=match.call())
  class(final.res) <- "SIGA"
  return(final.res)
}