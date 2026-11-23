# make all possible combinations of the variables to input the models
pop_models <- function(predictors=3,Rmax=c("obs", "est"), Rm=NULL, 
                       model=c("Verhulst","allee", "all"),
                       modelpars=c("yes","no"),data=data){
  
  
  predictors=3
  
  # required library
  require(stringr)
  require(stringi)
  require(minpack.lm)
  require(rsample)
  require(qpcR)
  require(plyr)
  
  
  allComb1 <- read.delim(file.path(PROJ, "data", "model_templates", "allComb2.txt"))
  allComb2 <- read.delim(file.path(PROJ, "data", "model_templates", "allCombAllee.txt"))
  
  # Create models from Verhulst (1838) in its Royama (2021) form
  
  # define model equations (uncompleted)
  if(Rmax=="est"){
    
    # Royama models
    Verhulst = "R~r-(c)*exp(Xt)"
    allee    = "R~r-(c)*exp(Xt)-(w)*exp(-Xt)"
    
    # Convert equations to Nt+1 = Nt
    Verhulst_N1a = c("r[j]-(c[j]", ")*exp(N1a)")
    Verhulst_N1b = c("r[j]-(c[j]", ")*exp(N1b)")
    
    allee_N1a = c("r[j]-(c[j]", ")*exp(N1a)-(w[j]", ")*exp(-N1a)")
    allee_N1b = c("r[j]-(c[j]", ")*exp(N1b)-(w[j]", ")*exp(-N1b)")
    
  }
  
  if(Rmax=="obs"){
    
    # Royama models
    Verhulst = str_glue("R~",max(Rm), "-(c)*exp(Xt)")
    allee    = str_glue("R~",max(Rm), "-(c)*exp(Xt)-(w)*exp(-Xt)")
    
    # Convert equations to Nt+1 = Nt
    Verhulst_N1a = c(str_glue(max(Rm),"-(c[j]"), ")*exp(N1a)")
    Verhulst_N1b = c(str_glue(max(Rm),"-(c[j]"), ")*exp(N1b)")
    
    allee_N1a = c(str_glue(max(Rm),"-(c[j]"), ")*exp(N1a)-(w[j]", ")*exp(-N1a)")
    allee_N1b = c(str_glue(max(Rm),"-(c[j]"), ")*exp(N1b)-(w[j]", ")*exp(-N1b)")
    
    # Verhulst = str_glue("R~",round(max(data$R),2), "-(c)*exp(Xt)")
    # allee    = str_glue("R~",round(max(data$R),2),"-(c)*exp(Xt)-(w)*exp(-Xt)")
    # 
    
    # Convert equations to Nt+1 = Nt
    Verhulst_N1a = c(str_glue(round(max(data$R),2),"-(c[j]"), ")*exp(N1a)")
    Verhulst_N1b = c(str_glue(round(max(data$R),2),"-(c[j]"), ")*exp(N1b)")
     
    allee_N1a = c(str_glue(round(max(data$R),2),"-(c[j]"), ")*exp(N1a)-(w[j]", ")*exp(-N1a)")
    allee_N1b = c(str_glue(round(max(data$R),2),"-(c[j]"), ")*exp(N1b)-(w[j]", ")*exp(-N1b)")
  }
  
  # add the k parameter to the Verhulst model
  intra <- unlist(stri_split_fixed(Verhulst, "c"),2)
  
  allee   <- stri_split_fixed(allee, c("c","w"),2)
  allee   <- unlist(stri_split_fixed(allee[[2]], "c",2))
  
  # loop to add variables to the model forms
  
  # for exp model
  expm <- c()
  for(i in 1:length(allComb1[,1])){
    expm[i]  <- str_glue(intra[1],"c",allComb1[[i,1]],intra[2],allComb1[[i,2]])
  }
  # for allee model
  expm2 <- c()
  for(i in 1:length(allComb2[,1])){
    expm2[i] <- str_glue(allee[1],"c",allComb2[[i,1]],allee[2],"w",allComb2[[i,2]],allee[3],allComb2[[i,3]])
  }
  
  #-------------------------------------------------------------------------------
  # load combinations for simulations
  compreds     <- read.delim(file.path(PROJ, "data", "model_templates", "allComb2_sim.txt"))
  compredsAlle <- read.delim(file.path(PROJ, "data", "model_templates", "allCombAllee_sim.txt"))
  
  # create the simulation formats for Royama normal model
  expmN1a <- c()
  expmN1b <- c()
  for(i in 1:length(compreds[,1])){
    expmN1a[i]  <- str_glue(Verhulst_N1a[1],compreds[[i,1]],Verhulst_N1a[2],compreds[[i,2]])
    expmN1b[i]  <- str_glue(Verhulst_N1b[1],compreds[[i,1]],Verhulst_N1b[2],compreds[[i,2]])
  }
  
  # create the simulation formats for Royama normal model
  expmN1a_alle <- c()
  expmN1b_alle <- c()
  for(i in 1:length(compredsAlle[,1])){
    expmN1a_alle[i]  <- str_glue(allee_N1a[1],compredsAlle[[i,1]],allee_N1a[2],compredsAlle[[i,2]],allee_N1a[3],compredsAlle[[i,3]])
    expmN1b_alle[i]  <- str_glue(allee_N1b[1],compredsAlle[[i,1]],allee_N1b[2],compredsAlle[[i,2]],allee_N1b[3],compredsAlle[[i,3]])
  }
  
  if(model=="Verhulst"){
    
    allmodels <- expm
    
    # simulation equation for Royama model
    simulabase     <- list(expmN1a,expmN1b)
  }
  
  if(model=="allee"){
    
    allmodels <- expm2
    
    # simulation equation for Allee model
    simulabase <- list(expmN1a_alle,expmN1b_alle)
  }
  
  # combine exp models and allee models in one array
  if(model=="all"){
    
    allmodels  <- c(expm,expm2)
    
    # simulation equation for all models
    simulabase <- list(c(expmN1a, expmN1a_alle), c(expmN1b, expmN1b_alle))
  }
  
  allComb3 <-  rbind(allComb1, allComb2)
  allComb3[is.na(allComb3)] <- ""
  #-----------------------------------------------------------------------------
  
  #extract the parameters and concatenate 
  allp  <- apply(allComb1, 1, paste, collapse="")
  allp2 <- apply(allComb2, 1, paste, collapse="")
  allp3 <- apply(allComb3, 1, paste, collapse="")
  
  # loop to extract only the letters of the parameters
  if(Rmax=="est"){
    rC  <- c("r","c") # add parameters for the Verhulst model 
    rCa <- c("r","c","w")
  }
  if(Rmax=="obs"){
    rC  <- c("c") # add parameters for the Verhulst model 
    rCa <- c("c","w")
  }
  
  allP  <- list()
  allP2 <- list()
  for(i in 1:length(allp)){
    allP[i]    <- str_extract_all(allp[i], '[[:lower:]]+')
    # append the base parameters
    allP2[[i]] <- append(rC,unlist(allP[i])) # for
  }
  
  # loop to extract only the letters of the parameters
   # add parameters for the Verhulst-allee model
  allP  <- list()
  allP3 <- list()
  for(i in 1:length(allp2)){
    allP[i]    <- str_extract_all(allp2[i], '[[:lower:]]+')
    # append the base parameters
    allP3[[i]] <- append(rCa,unlist(allP[i])) 
  }
  
  # combine all parameters for the Verhulst-allee model
  allpar <- c(allP2, allP3)
  
  #-----------------------------------------------------------------------------
  # insert "=" and numbers to generate starting values
  
  # loop to add comas and initial value to the parameters for the Verhulst model 
  ppp  <- list()
  pars1 <- list() 
  for(i in 1:length(allP2)){
    x = c(unlist(allP2[[i]])) 
    pp <- paste(x,1,sep = '=',collapse = ' , ')
    ppp[[i]]  <- paste(pp,sep="")
    pars1[[i]] <- eval(parse(text = paste("list(",ppp[[i]],")")))
  }
  
  # loop to add comas and initial value to the parameters for the Verhulst.allee 
  ppp  <- list()
  pars2 <- list() 
  for(i in 1:length(allP3)){
    x = c(unlist(allP3[[i]])) 
    pp <- paste(x,1,sep = '=',collapse = ' , ')
    ppp[[i]]  <- paste(pp,sep="")
    pars2[[i]] <- eval(parse(text = paste("list(",ppp[[i]],")")))
  }
  
  # loop to add comas and initial value to the parameters for the Verhulst + allee 
  ppp  <- list()
  pars3 <- list() 
  for(i in 1:length(allpar)){
    x = c(unlist(allpar[[i]])) 
    pp <- paste(x,1,sep = '=',collapse = ' , ')
    ppp[[i]]  <- paste(pp,sep="")
    pars3[[i]] <- eval(parse(text = paste("list(",ppp[[i]],")")))
  }
  #-----------------------------------------------------------------------------
  # make list of models and parameters for return
  
  # for Verhulst
  ver_mods_pars      <- list(expm, pars1)
  
  # for Verhulst-allee
  alle_mods_pars     <- list(expm2, pars2)
  
  # for Verhulst-allee
  Ver_alle_mods_pars <- list(allmodels, pars3)
  
  # return lists: 1. models; 2. parameters
  # if(model=="Verhulst"){return(ver_mods_pars)}
  # if(model=="allee"){return(alle_mods_pars)}
  # 
  if(model=="Verhulst"){allmodels <- ver_mods_pars}
  if(model=="allee")   {allmodels <- alle_mods_pars}
  if(model=="all")     {allmodels <- Ver_alle_mods_pars}
  #-------------------------------------------------------------------------------
  # run all models combinations
  models <- list()
  for(i in 1:length(allmodels[[1]])){
    models[[i]]  <- nlsLM(allmodels[[1]][i], 
                          start = allmodels[[2]][[i]], 
                          control = nls.lm.control(maxiter = 1000), 
                          algorithm = "BFGS",
                          data = data)
  }
  
  # function to calculate R2
  R2 <- function(modelo){
    (cor(data$R,predict(modelo))^2)
  }
  
  # get all models in a list for final table creation and selection
  d.2     <- sapply(models, R2)
  # AICc from qpcR: k = number of estimated coefficients (MuMIn::AICc counts sigma as well)
  all.AIC <- sapply(models, qpcR::AICc)
  all.BIC <- sapply(models, BIC)
  forms   <- sapply(models, formula)
  params  <- sapply(models, coef)
  loglike <- sapply(models, logLik)
  npars   <- sapply(params, length)
  w       <- akaike.weights(all.AIC)
  
  # loop to extract parameters values
  par.df  <- list()
  for(i in 1:length(params)){
    par.df[[i]] <- as.data.frame(t(params[[i]]))
  }
  
  # create data frame with all model results
  all.models <- data.frame(Models=as.character(forms),
                           N.param=npars,
                           AICc=all.AIC, 
                           DAICc=w$deltaAIC,
                           # BIC=all.BIC,
                           wi=w$weights,
                           R2=d.2)
  
  # create data frame with all model results with parameters
  all.models1 <- data.frame(Models=as.character(forms),
                            N.param=npars,
                            rbind.fill(par.df), # add parameter values
                            AICc=all.AIC,
                            DAICc=w$deltaAIC,
                            # BIC=all.BIC,
                            wi=w$weights,
                            R2=d.2)
  
  all.models  <- list(all.models, models, simulabase)
  all.models1 <- list(all.models1,models, simulabase)
  
  if(modelpars=="no") {return(all.models)}
  if(modelpars=="yes"){return(all.models1)}
  
}
