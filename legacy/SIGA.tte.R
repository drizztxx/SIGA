SIGA.tte <- function(data,
                        trt="treat=1",
                        ctrl="treat=2",
                        status="status",
                        time="time",
                        varScreen=c(1,5),
                        varused=NULL,
                        popSize=100,
                        objFunction=c("objF.cox"),
                        MAXGEN = 400,        ## max Number of generations
                        GGAP = 0.9,          ## Generation gap, how many new individuals are created
                        SEL_F = 'sus',       ## Name of selection function
                        XOV_F = 'reclin',     ## Name of recombination function for individuals
                        MUT_F = 'mutbga',       ## Name of mutation function for individuals
                        cutoff = 10,
                        verbose = TRUE,
                        ...
){
  library(randomForestSRC)
  warnold <- getOption("warn")
  options(warn=-1)
  trtvar <- strsplit(trt,'=')[[1]][1]
  if (!is.null(varused)){
    col.vim <- varused
    top <- NROW(col.vim)
  } else if (varScreen[1]) { ## Screen the top important variables via survival RF
    top <- varScreen[2]
    screen.res <- VariableScreen(data=data,trt=trt,ctrl=ctrl,status = status,time=time,varScreen = varScreen,resptyp = "tte")
    col.vim <- screen.res$col.vim[1:top]
  } else {
    top <- NCOL(data) - 3
    col.vim <- colnames(data)[!colnames(data) %in% c(trtvar,time,status)]
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
  Chrom <- crtrp(popSize,FieldDR)
  Chrom.pool <- list()
  Chrom.pool[[1]] <- Chrom
  
  ## Decoding the generated Chrom according to the genetic representation rule.
  if (top == 1){
    decoddata <- as.data.frame(data[,col.vim])
    colnames(decoddata) <- col.vim
  } else decoddata <- data[,col.vim]
  decod.subset <- decodsubset(data=decoddata,Chrom = Chrom)
  subset.pool <- list()
  subset.pool[[1]] <- decod.subset$subgroup.desc
  
  ## Caculate the treatment effect (that is the pvalue caculated through the glm model) for every subgroups
  ObjV <- do.call(objFunction,list(subgroups=decod.subset$subgroup.desc, 
                                   data=data[,c(col.vim,time,status,trtvar)],trt = trt,
                                   ctrl = ctrl,status=status,time=time,cutoff=cutoff,... ))
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
    FitnV = ranking(objv);
    
    ## Select individuals from population
    SelCh = gatbxr::select(SEL_F, Chrom, FitnV, GGAP);

    ## Recombine selected individuals (crossover)
    SelCh=reclin(SelCh);
    
    ## Mutate offspring
    SelCh=gatbxr::mutate(MUT_F, SelCh, FieldDR = FieldDR);
    
    ##Evaluate offspring, call objective function
    ObjVSel <- do.call(objFunction,list(subgroups=decodsubset(data=data[,col.vim],Chrom = SelCh)$subgroup.desc, 
                                        data=data[,c(col.vim,time,status,trtvar)],trt = trt,
                                        ctrl = ctrl,status=status,time=time,cutoff=cutoff,... ))
    
    
    ## Insert offspring in population replacing parents
    rs = reins(Chrom, SelCh,1,c(1,1),objv,ObjVSel$objv);
    Chrom = rs$Chrom;
    objv  = rs$ObjVCh
    
    gen = gen + 1;
    Chrom.pool[[gen]] <- Chrom
    objv.pool[[gen]] <- objv
    subset.pool[[gen]] <- decodsubset(data=decoddata,Chrom = Chrom)$subgroup.desc
  }
  ## End of script
  options(warn=warnold)
  dsn <- deparse(substitute(data))
  if (!exists("screen.res")) screen.res <- NULL
  final.res <- list(Chrom=Chrom.pool,Objv=objv.pool,Subset=subset.pool,Best=Best,varused=col.vim,dsn=dsn,resptyp="tte",
                    screen.res=screen.res,call=match.call())
  class(final.res) <- "SIGA"
  return(final.res)
}