#################################################################################
############# Data that one can mock on gwsemPLS ################################
#################################################################################
gwsemPLS_demoDataLink<-"https://disk.pku.edu.cn/link/AA4895174AFA764ABF8E1B3864BDFFE3DC"


#################################################################################
############# Install and load packages needed   ################################
#################################################################################
# Function to check and install/load packages
install_and_load_package <- function(package_name) {
  # Check if the package is installed
  if (!requireNamespace(package_name, quietly = TRUE)) {
    # Install the package if not installed
    install.packages(package_name)
  }
  # Load the package
  library(package_name, character.only = TRUE)
}


##################################################################################################################################################################
#Part 1. Basic Unit of the gwsemPLS package 
#When your dataset contains only ONE snp and other variables/indicators, you can examine the model fit, reliability, coefficient estimates, and their stanardard errors.
#Note that standard errors and confidence intervals are obtained via Bootstrapping and your should specify the "bootstrap=T" to get the corresponding results.
#The results in the format of R objects of the functions below are termed as "SingleSnpObj".
#The "SingleSnpObj" is a list where SingleSnpObj[[1]] is the model fit results without bootstrapping and SingleSnpObj[[2]] is the bootstrapped results.
#The "SingleSnp_summary(SingleSnpObj)" can yield summary statistics and information about the model specified in fitting gwsem.
#The "gwsemPLS_plot(SingleSnpObj)" can show structure of the model specified in fitting gwsem.
#(NOT RUN)
#SingleSnpObj<-SingleSnp_gwSEM_pls_StadGWAS(phenoData_with_singleSnp, depVar = c("tobacco"), covariates=c('pc1','pc2','pc3','pc4','pc5'))
#SingleSnpObj<-SingleSnp_gwSEM_pls_SingleFact(phenoData_with_singleSnp,itemNames = c("tobacco", "cannabis", "alcohol"),covariates = c("pc1", "pc2", "pc3", "pc4", "pc5"))
#SingleSnpObj<-SingleSnp_gwSEM_pls_GxE(phenoData_with_singleSnp,depVar  = "phe",covariates = c("mod", "pc1", "pc2", "pc3", "pc4","pc5"), gxe = "mod")
#SingleSnpObj<-SingleSnp_gwSEM_pls_SingleFactRes(phenoData_with_singleSnp, itemNames = c("tobacco", "cannabis", "alcohol"), res = c("tobacco", "cannabis", "alcohol"), covariates=c('pc1','pc2','pc3','pc4','pc5'))
#SingleSnpObj<-SingleSnp_gwSEM_pls_MultiFact(phenoData_with_singleSnp,FitemNames=list(c("A1", "A2", "A3"), c("B1", "B2", "B3")),covariates= c('pc1','pc2','pc3','pc4','pc5')) 
#SingleSnp_summary(SingleSnpObj)
##################################################################################################################################################################
#Part2. Multiple SNPs are available
#When your dataset contains multiple snps and other variables, you can obtain the results of parameter estimates and their pvalues
#The same model specified is applied to all SNPs and other variables/indicators
#The results in the format of R objects of the functions below are termed as "MultipleSnpObj".
#The "MultipleSnpObj" is a dataframe containing (path) coefficient estimates and their pvalues in the form of "snpname: pathHead -> pathTail"
#The "detect_snp_rows(MultipleSnpObj)" can filter non-snp effects and deliver the SNP and table results containing "SNP" and "P" 
#(NOT RUN)
#MultipleSnpObj<-gwSEM_pls_StadGWAS(phenoData_with_manySnp, depVar = c("tobacco"), covariates=c('pc1','pc2','pc3','pc4','pc5'))
#MultipleSnpObj<-gwSEM_pls_SingleFact(phenoData_with_manySnp,itemNames = c("tobacco", "cannabis", "alcohol"),covariates = c("pc1", "pc2", "pc3", "pc4", "pc5"))
#MultipleSnpObj<-gwSEM_pls_GxE(phenoData_with_manySnp,depVar  = "phe",covariates = c("mod", "pc1", "pc2", "pc3", "pc4","pc5"), gxe = "mod")
#MultipleSnpObj<-gwSEM_pls_SingleFactRes(phenoData_with_manySnp, itemNames = c("tobacco", "cannabis", "alcohol"), res = c("tobacco", "cannabis", "alcohol"), covariates=c('pc1','pc2','pc3','pc4','pc5'))
#MultipleSnpObj<-gwSEM_pls_MultiFact(phenoData_with_manySnp,FitemNames=list(c("A1", "A2", "A3"), c("B1", "B2", "B3")),covariates= c('pc1','pc2','pc3','pc4','pc5')) 
#snpEffectOnly_MultipleSnpObj<-detect_snp_rows(MultipleSnpObj)
##################################################################################################################################################################
#Part3. Visualize the pvalue as Manhattan plots
#Use the "snpEffectOnly_MultipleSnpObj" as input to execute "gwsemPLS_manhattanPlot()"
#If chromosome and position information is also available, you can add two columns- "CHR" and ”BP“- to "snpEffectOnly_MultipleSnpObj".
#The "suggestiveline" and "genomewideline" configurations can be changed in the "gwsemPLS_manhattanPlot()" 
#(NOT RUN)
#gwsemPLS_manhattanPlot(snpEffectOnly_MultipleSnpObj)
##################################################################################################################################################################



#################################################################################
#############     Supportive Functions       ####################################
#################################################################################

# Make sure the dataset has SNP/snp in column names
checkSnpVariable<-function(dataset,returnDat='F'){
  # (1) Change all column names to lower case
  colnames(dataset) <- tolower(colnames(dataset))
  # (2) Check if there is at least one column name containing 'snp'
  snp_columns <- grep("snp", colnames(dataset), value = TRUE)
  if (length(snp_columns) == 0) {
    print("No SNP variables are found!")
    return(NULL)
  }
  if(returnDat){
    dataset
  }
}

# Identify column names' indexes of snp and non-snp variables
processForIteration <- function(dataset) {
  dataset<-checkSnpVariable(dataset,T) 
  #Identify columns named with 'snp' as a substring
  snp_identifier <- grep("snp", colnames(dataset))
  #Identify columns that do not contain 'snp'
  non_snp_identifier  <- setdiff(1:ncol(dataset), snp_identifier)
  # Return the modified dataset
  process_data<-list()
  process_data[[1]]<-snp_identifier
  process_data[[2]]<-non_snp_identifier
  process_data
}

# Custom Manhattan plot function
custom_simplified_manhattan <- function(gwasResults, suggestiveline, genomewideline,...) {
  # Ensure the required columns are present
  if (!all(c('SNP', 'P') %in% colnames(gwasResults))) {
    stop("Missing required columns in gwasResults")
  }
  # Extract p-values
  p_values <- gwasResults[['P']]
  # Calculate -log10(p-values)
  log_p_values <- -log10(p_values)
  # Create a sequence for the x-axis (SNP positions)
  snp_positions <- seq_along(gwasResults[['SNP']])
  # Determine colors based on the suggestiveline
  point_colors <- ifelse(log_p_values > suggestiveline, "red", "green")
  # Create the plot
  plot(snp_positions, log_p_values, 
       type = "p", 
       pch = 20, 
       col = point_colors, 
       xlab = "SNP Index", 
       ylab = "-log10(P-value)", 
       main = "Manhattan Plot", 
       ...)
  abline(h = suggestiveline, col = "black", lty = 2)
}
gwsemPLS_manhattanPlot<-function(gwasResults,suggestiveline = -log10(0.05), genomewideline = -log10(5e-8)){
  print('Showing chromosome and position information is permissible:')
  print('Add rows named \'CHR\' and \'BP\' to the dataframe of the manhattanPlot()\'s input.')
  if (ncol(gwasResults)>2){
    (manhattan(gwasResults, suggestiveline=suggestiveline, genomewideline=genomewideline))
  }else{
    (custom_simplified_manhattan(gwasResults, suggestiveline, genomewideline))
  }
}

# Detect Pvalue results and extract rows with two SNP/snp at their names 
# Define the function
detect_snp_rows <- function(data) {
  snpRownames<-rownames(data)
  # Use sapply to apply a function to each element of the vector
  # Define the function
  T_F.selector<-sapply(snpRownames, function(row) {
    # Count the occurrences of 'snp' (case insensitive) in the entire string
    count <- str_count(tolower(row), "snp")
    # Return TRUE if 'snp' appears at least twice
    count >= 2
  })  
  res<-data[T_F.selector,]
  res<-cbind(rownames(res),res$pvalues)
  colnames(res)<-c('SNP','P')
  
  res<-as.data.frame(res)
  res$P<-as.numeric(res$P)
  res
  #gwasResults<-detect_snp_rows((gwSEM_pls_StadGWAS_ResJbject))
}

# Extract pvalue from bootstrapped seminr object
extractPvalue<-function(bootedObj){ 
  #extractPvalue(bootedObj)
  boot_summary<-summary(bootedObj)
  # gather paths and t-values
  paths <- boot_summary$bootstrapped_paths[, "Original Est."]
  tvalues <- boot_summary$bootstrapped_paths[, "T Stat."]
  # degrees of freedom will be the number of rows in the data sample - number of parameters
  df = nrow(bootedObj$data) 
  # calculate pvalues from tvalues and df; round to 3 decimal places. multiplied by 2 to get the two-tailed p-value
  pvalues <- round( 2*pt(tvalues, df, lower.tail = FALSE), 3)
  # make a table of paths, tvalues, pvalues
  data.frame(paths, tvalues, pvalues)
}

# Get summary resutls for one 
SingleSnp_summary<-function(SingleSnp_gwSEM_plsObj){
  summary(SingleSnp_gwSEM_plsObj[[1]])
}

# Plot the model
gwsemPLS_plot<-function(SingleSnp_gwSEM_plsObj){
  if(length(SingleSnp_gwSEM_plsObj)>1){
    plot(SingleSnp_gwSEM_plsObj[[2]], title = "Bootstrapped Model")
  }else{
    plot(SingleSnp_gwSEM_plsObj[[1]], title = "Model Without Bootstrapping")
  }
}


#################################################################################
######################  Standard GWAS ###########################################
#################################################################################

SingleSnp_gwSEM_pls_StadGWAS <- function(dat, depVar, covariates,bootstrap=F,nboot = 200,cores = 2,seed = 123) {
  # List of packages to check and install/load
  packages_to_check <- c("seminr", "OpenMx", "gwsem","qqman","stringr","plotly")
  # Loop through the packages and install/load them
  for (pkg in packages_to_check) {
    install_and_load_package(pkg)
  }
  checkSnpVariable(dat) 
  # Ensure depVar is treated as a single item, not a list
  if (length(depVar) > 1) {
    stop("depVar should contain only one dependent variable.")
  }
  
  # Create the composite for the dependent variable
  depVarComposite <- paste0("composite('", toupper(depVar), "',single_item('", depVar, "'))")
  
  # Create the SNP composite
  snpComposite <- "composite('SNP',single_item('snp'))"
  
  # Create composites for covariates
  covariateComposites <- sapply(covariates, function(cov) {
    paste0("composite('", toupper(cov), "',single_item('", cov, "'))")
  }, USE.NAMES = FALSE)
  
  # Combine all measurements
  measurements <- paste("measurements <- constructs(",
                        depVarComposite, ",",
                        snpComposite, ",",
                        paste(covariateComposites, collapse = ",
                         "),
                        ")")
  # Create structural relationships for the dependent variable with SNP and each covariate
  structuralPaths <- sapply(covariates, function(cov) {
    paste0("paths(from = '", toupper(cov), "', to = c('", toupper(depVar), "'))")
  }, USE.NAMES = FALSE)
  
  # Add the path from SNP to the dependent variable
  structuralPaths <- c(paste0("paths(from = 'SNP', to = c('", toupper(depVar), "'))"), structuralPaths)
  
  # Combine the structural relationships
  structural <- paste("structural <- relationships(",
                      paste(structuralPaths, collapse = ","),")")
  
  
  eval(parse(text=measurements))
  eval(parse(text=structural))
  beforeBootRes<-estimate_pls(data = dat,
               measurement_model = measurements,
               structural_model = structural,
               inner_weights = path_weighting)
  finalRes<-list()
  finalRes[[1]]<-beforeBootRes
  if(bootstrap){
    finalRes[[2]]<-bootstrap_model(beforeBootRes,nboot = nboot,cores = cores,seed=seed)
  }
  finalRes
}
gwSEM_pls_StadGWAS <- function(dat, depVar, covariates, nboot = 200,cores = 2,seed = 123) {
  #Before Iteration
  snpColID<-processForIteration(dat)[[1]]
  snpNames<-names(dat)[snpColID]
  nonsnpColID<-processForIteration(dat)[[2]]
  
  pathPvalues<-NULL
  #Start Iteration for each snp
  for (i in 1:length(snpColID)){
    datForIteration<-cbind(dat[,snpColID[i]],dat[,nonsnpColID])
    names(datForIteration)[1]<-'snp'
    SingleSnpRes<-SingleSnp_gwSEM_pls_StadGWAS(datForIteration, depVar = depVar, covariates=covariates)
    
    
    boot_mobi_pls<-bootstrap_model(SingleSnpRes[[1]],nboot = nboot,cores = cores,seed=seed)
    temp_pathPvalues<-extractPvalue(boot_mobi_pls)
    temp_snpName<-snpNames[i]
    rownames(temp_pathPvalues)<-paste(temp_snpName,rownames(temp_pathPvalues),sep=': ')
    pathPvalues<-rbind(pathPvalues,temp_pathPvalues)
  }
  pathPvalues
}

#################################################################################
##################### One factor model ##########################################
#################################################################################

SingleSnp_gwSEM_pls_SingleFact <- function(dat,itemNames, covariates,bootstrap=F,nboot = 200,cores = 2,seed = 123) {
  # List of packages to check and install/load
  packages_to_check <- c("seminr", "OpenMx", "gwsem","qqman","stringr","plotly")
  # Loop through the packages and install/load them
  for (pkg in packages_to_check) {
    install_and_load_package(pkg)
  }
  checkSnpVariable(dat)
  # Create the measurements
  measurements <- paste0(
    "measurements <- constructs(",
    "reflective('F',c('", paste(itemNames, collapse = "','"), "')),",
    paste0(
      "composite('", toupper(itemNames), "',single_item('", itemNames, "'))",
      collapse = ","
    ),
    ",",
    "composite('SNP',single_item('snp')),",
    paste0(
      "composite('", toupper(covariates), "',single_item('", covariates, "'))",
      collapse = ","
    ),
    ")"
  )
  
  # Create the structural relationships
  structural <- paste0(
    "structural <- relationships(",
    "paths(from = 'SNP', to = c('F')),",
    paste0(
      "paths(from = '", toupper(covariates), "', to = c('", toupper(itemNames), "'))",
      collapse = ","
    ),
    ")"
  )
  
  eval(parse(text=measurements))
  eval(parse(text=structural))
  beforeBootRes<-estimate_pls(data = dat,
                              measurement_model = measurements,
                              structural_model = structural,
                              inner_weights = path_weighting)
  finalRes<-list()
  finalRes[[1]]<-beforeBootRes
  if(bootstrap){
    finalRes[[2]]<-bootstrap_model(beforeBootRes,nboot = nboot,cores = cores,seed=seed)
  }
  finalRes
}
gwSEM_pls_SingleFact <- function(dat, itemNames, covariates, nboot = 200,cores = 2,seed = 123) {
  #Before Iteration
  snpColID<-processForIteration(dat)[[1]]
  snpNames<-names(dat)[snpColID]
  nonsnpColID<-processForIteration(dat)[[2]]
  pathPvalues<-NULL
  #Start Iteration for each snp
  for (i in 1:length(snpColID)){
    datForIteration<-cbind(dat[,snpColID[i]],dat[,nonsnpColID])
    names(datForIteration)[1]<-'snp'
    SingleSnpRes<-SingleSnp_gwSEM_pls_SingleFact(datForIteration, itemNames = itemNames, covariates=covariates)
    boot_mobi_pls<-bootstrap_model(SingleSnpRes[[1]],nboot = nboot,cores = cores,seed=seed)
    temp_pathPvalues<-extractPvalue(boot_mobi_pls)
    temp_snpName<-snpNames[i]
    rownames(temp_pathPvalues)<-paste(temp_snpName,rownames(temp_pathPvalues),sep=': ')
    pathPvalues<-rbind(pathPvalues,temp_pathPvalues)
  }
  pathPvalues
}

#################################################################################
############  Gene environment interaction  #####################################
#################################################################################

SingleSnp_gwSEM_pls_GxE <- function(dat, depVar, covariates, gxe,bootstrap=F,nboot = 200,cores = 2,seed = 123) {
  # List of packages to check and install/load
  packages_to_check <- c("seminr", "OpenMx", "gwsem","qqman","stringr","plotly")
  # Loop through the packages and install/load them
  for (pkg in packages_to_check) {
    install_and_load_package(pkg)
  }
  checkSnpVariable(dat)
  # Create the composite for the dependent variable and gxe
  depVarComposite <- paste0("composite('", toupper(depVar), "',single_item('", depVar, "'))")
  gxeComposite <- paste0("composite('", toupper(gxe), "',single_item('", gxe, "'))")
  
  # Create composites for covariates
  covariateComposites <- paste0("composite('", toupper(covariates), "',single_item('", covariates, "'))", collapse = ",")
  
  # Create the SNP composite
  snpComposite <- "composite('SNP',single_item('snp'))"
  
  # Create the interaction term
  interactionTerm <- paste0("interaction_term(iv = 'SNP', moderator = '", toupper(gxe), "', method = product_indicator)")
  
  # Combine all measurements
  measurements <- paste("measurements <- constructs(",
                        gxeComposite, ",",
                        depVarComposite, ",",
                        snpComposite, ",",
                        covariateComposites, ",",
                        interactionTerm,
                        ")")
  
  # Create structural relationships
  structuralPaths <- paste0("paths(from = '", toupper(covariates), "', to = c('", toupper(depVar), "'))", collapse = ",")
  structural <- paste("structural <- relationships(",
                      "paths(from = 'SNP', to = c('", toupper(depVar), "')),",
                      "paths(from = '", toupper(gxe), "', to = c('", toupper(depVar), "')),",
                      structuralPaths, ",",
                      "paths(from = 'SNP*", toupper(gxe), "', to = '", toupper(depVar), "')",
                      ")",sep='')
  
  
  
  
  eval(parse(text=measurements))
  eval(parse(text=structural))
  
  # Define the function to remove duplicate rows from a matrix
  removeDuplicateRows <- function(matrix) {
    # Convert the matrix to a data frame for easy handling
    df <- as.data.frame(matrix)
    
    # Remove duplicate rows
    unique_df <- unique(df)
    
    # Convert back to a matrix, if needed
    unique_matrix <- as.matrix(unique_df)
    
    # Return the unique matrix
    return(unique_matrix)
  }
  structural<-removeDuplicateRows(structural)
  beforeBootRes<-estimate_pls(data = dat,
                              measurement_model = measurements,
                              structural_model = structural,
                              inner_weights = path_weighting)
  finalRes<-list()
  finalRes[[1]]<-beforeBootRes
  if(bootstrap){
    finalRes[[2]]<-bootstrap_model(beforeBootRes,nboot = nboot,cores = cores,seed=seed)
  }
  finalRes
}
gwSEM_pls_GxE <- function(dat, depVar,gxe, covariates, nboot = 200,cores = 2,seed = 123) {
  #Before Iteration
  snpColID<-processForIteration(dat)[[1]]
  snpNames<-names(dat)[snpColID]
  nonsnpColID<-processForIteration(dat)[[2]]
  pathPvalues<-NULL
  #Start Iteration for each snp
  for (i in 1:length(snpColID)){
    datForIteration<-cbind(dat[,snpColID[i]],dat[,nonsnpColID])
    names(datForIteration)[1]<-'snp'
    SingleSnpRes<-SingleSnp_gwSEM_pls_GxE(datForIteration, depVar = depVar, covariates=covariates,gxe=gxe)
    boot_mobi_pls<-bootstrap_model(SingleSnpRes[[1]],nboot = nboot,cores = cores,seed=seed)
    temp_pathPvalues<-extractPvalue(boot_mobi_pls)
    temp_snpName<-snpNames[i]
    rownames(temp_pathPvalues)<-paste(temp_snpName,rownames(temp_pathPvalues),sep=': ')
    pathPvalues<-rbind(pathPvalues,temp_pathPvalues)
  }
  pathPvalues
}

#################################################################################
##############   One factor residuals model #####################################
#################################################################################

SingleSnp_gwSEM_pls_SingleFactRes <- function(dat,itemNames, res, covariates,bootstrap=F,nboot = 200,cores = 2,seed = 123) {
  # List of packages to check and install/load
  packages_to_check <- c("seminr", "OpenMx", "gwsem","qqman","stringr","plotly")
  # Loop through the packages and install/load them
  for (pkg in packages_to_check) {
    install_and_load_package(pkg)
  }
  checkSnpVariable(dat)
  # Create the reflective construct for the factor with specified items
  reflectiveConstruct <- paste0("reflective('F',c('", paste(itemNames, collapse = "','"), "'))")
  
  # Create composites for each item in res
  itemComposites <- sapply(res, function(item) {
    paste0("composite('", toupper(item), "',single_item('", item, "'))")
  }, USE.NAMES = FALSE)
  
  # Create the SNP composite
  snpComposite <- "composite('SNP',single_item('snp'))"
  
  # Create composites for covariates
  covariateComposites <- sapply(covariates, function(cov) {
    paste0("composite('", toupper(cov), "',single_item('", cov, "'))")
  }, USE.NAMES = FALSE)
  
  # Combine all measurements
  measurements <- paste("measurements <- constructs(",
                        reflectiveConstruct, ",",
                        paste(c(itemComposites, snpComposite, covariateComposites), collapse = ","),")")
  
  # Create structural relationships for each item and covariate
  structuralItems <- sapply(itemNames, function(item) {
    paste0("paths(from = 'SNP', to = c('", toupper(item), "')),",
           paste0("paths(from = '", toupper(covariates), "', to = c('", toupper(item), "'))", collapse = ","))
  }, USE.NAMES = FALSE)
  
  # Combine the structural relationships
  structural <- paste("structural <- relationships(",
                      paste(structuralItems, collapse = ","),
                      ")")
  
  eval(parse(text=measurements))
  eval(parse(text=structural))
  
  # Define the function to remove duplicate rows from a matrix
  removeDuplicateRows <- function(matrix) {
    # Convert the matrix to a data frame for easy handling
    df <- as.data.frame(matrix)
    
    # Remove duplicate rows
    unique_df <- unique(df)
    
    # Convert back to a matrix, if needed
    unique_matrix <- as.matrix(unique_df)
    
    # Return the unique matrix
    return(unique_matrix)
  }
  structural<-removeDuplicateRows(structural)
  beforeBootRes<-estimate_pls(data = dat,
                              measurement_model = measurements,
                              structural_model = structural,
                              inner_weights = path_weighting)
  finalRes<-list()
  finalRes[[1]]<-beforeBootRes
  if(bootstrap){
    finalRes[[2]]<-bootstrap_model(beforeBootRes,nboot = nboot,cores = cores,seed=seed)
  }
  finalRes
}
gwSEM_pls_SingleFactRes <- function(dat, itemNames,res, covariates, nboot = 200,cores = 2,seed = 123) {
  #Before Iteration
  snpColID<-processForIteration(dat)[[1]]
  snpNames<-names(dat)[snpColID]
  nonsnpColID<-processForIteration(dat)[[2]]
  pathPvalues<-NULL
  #Start Iteration for each snp
  for (i in 1:length(snpColID)){
    datForIteration<-cbind(dat[,snpColID[i]],dat[,nonsnpColID])
    names(datForIteration)[1]<-'snp'
    SingleSnpRes<-SingleSnp_gwSEM_pls_SingleFactRes(datForIteration, itemNames = itemNames, res=res,covariates=covariates)
    boot_mobi_pls<-bootstrap_model(SingleSnpRes[[1]],nboot = nboot,cores = cores,seed=seed)
    temp_pathPvalues<-extractPvalue(boot_mobi_pls)
    temp_snpName<-snpNames[i]
    rownames(temp_pathPvalues)<-paste(temp_snpName,rownames(temp_pathPvalues),sep=': ')
    pathPvalues<-rbind(pathPvalues,temp_pathPvalues)
  }
  pathPvalues
}

#################################################################################
########################## Multi-factor model ###################################
#################################################################################

SingleSnp_gwSEM_pls_MultiFact <- function(dat,FitemNames, covariates,bootstrap=F,nboot = 200,cores = 2,seed = 123) {
  # List of packages to check and install/load
  packages_to_check <- c("seminr", "OpenMx", "gwsem","qqman","stringr","plotly")
  # Loop through the packages and install/load them
  for (pkg in packages_to_check) {
    install_and_load_package(pkg)
  }
  checkSnpVariable(dat)
  # Validate FitemNames is a list
  if (!is.list(FitemNames)) {
    stop("FitemNames must be a list.")
  }
  
  # Number of factors
  numFactors <- length(FitemNames)
  
  # Create reflective constructs for each factor
  reflectiveConstructs <- lapply(seq_len(numFactors), function(idx) {
    paste0("reflective('F", idx, "',c('", paste(FitemNames[[idx]], collapse = "','"), "'))")
  })
  
  # Create composites for each item
  itemComposites <- unique(unlist(lapply(FitemNames, function(items) {
    sapply(items, function(item) {
      paste0("composite('", toupper(item), "',single_item('", item, "'))")
    })
  })))
  
  # Create the SNP composite
  snpComposite <- "composite('SNP',single_item('snp'))"
  
  # Create composites for covariates
  covariateComposites <- sapply(covariates, function(cov) {
    paste0("composite('", toupper(cov), "',single_item('", cov, "'))")
  }, USE.NAMES = FALSE)
  
  # Combine all measurements
  measurements <- paste("measurements <- constructs(",
                        paste(c(reflectiveConstructs, itemComposites, snpComposite, covariateComposites), collapse = ","),")")
  
  # Structural relationships
  # Paths from SNP to factors
  snpToFactorsPaths <- paste0("paths(from = 'SNP', to = c(", paste(sprintf("'F%d'", 1:numFactors), collapse = ","), "))")
  
  # Paths from covariates to all items
  allItems <- unique(unlist(FitemNames))
  covariatePaths <- sapply(covariates, function(cov) {
    paste0("paths(from = '", toupper(cov), "', to = c('", paste(toupper(allItems), collapse = "','"), "'))")
  }, USE.NAMES = FALSE)
  
  # Combine the structural relationships
  structural <- paste("structural <- relationships(",
                      snpToFactorsPaths, ",",
                      paste(covariatePaths, collapse = ","),")")
  
  eval(parse(text=measurements))
  eval(parse(text=structural))
  
  # Define the function to remove duplicate rows from a matrix
  removeDuplicateRows <- function(matrix) {
    # Convert the matrix to a data frame for easy handling
    df <- as.data.frame(matrix)
    
    # Remove duplicate rows
    unique_df <- unique(df)
    
    # Convert back to a matrix, if needed
    unique_matrix <- as.matrix(unique_df)
    
    # Return the unique matrix
    return(unique_matrix)
  }
  structural<-removeDuplicateRows(structural)
  beforeBootRes<-estimate_pls(data = dat,
                              measurement_model = measurements,
                              structural_model = structural,
                              inner_weights = path_weighting)
  finalRes<-list()
  finalRes[[1]]<-beforeBootRes
  if(bootstrap){
    finalRes[[2]]<-bootstrap_model(beforeBootRes,nboot = nboot,cores = cores,seed=seed)
  }
  finalRes
}
gwSEM_pls_MultiFact <- function(dat, FitemNames, covariates, nboot = 200,cores = 2,seed = 123) {
  #Before Iteration
  snpColID<-processForIteration(dat)[[1]]
  snpNames<-names(dat)[snpColID]
  nonsnpColID<-processForIteration(dat)[[2]]
  pathPvalues<-NULL
  #Start Iteration for each snp
  for (i in 1:length(snpColID)){
    datForIteration<-cbind(dat[,snpColID[i]],dat[,nonsnpColID])
    names(datForIteration)[1]<-'snp'
    SingleSnpRes<-SingleSnp_gwSEM_pls_MultiFact(datForIteration, FitemNames = FitemNames, covariates=covariates)
    boot_mobi_pls<-bootstrap_model(SingleSnpRes[[1]],nboot = nboot,cores = cores,seed=seed)
    temp_pathPvalues<-extractPvalue(boot_mobi_pls)
    temp_snpName<-snpNames[i]
    rownames(temp_pathPvalues)<-paste(temp_snpName,rownames(temp_pathPvalues),sep=': ')
    pathPvalues<-rbind(pathPvalues,temp_pathPvalues)
  }
  pathPvalues
}
