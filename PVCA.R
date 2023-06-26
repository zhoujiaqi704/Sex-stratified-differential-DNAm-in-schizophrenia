#' Principal Variance Component Analysis
#'
#' This function perform PVCA without interaction.
#'
#' @param data Data matrix
#' @param trait Sample information used in the liner model
#' @param pct_threshold The threshold for PC Defaults to .5876
#' @param main The title in plot Defaults to "PVCA plot"
#' @param xlab The xlab in plot Defaults to "Effects"
#' @keywords PVCA
#' @return PVCA plot
#' @export
#' @examples
#' pvca(data, trait)
#' pvca2(data, trait)
#' @export

pvca <- function(data, trait, pct_threshold = .5876, main="PVCA plot", xlab="Effects"){
    pct_threshold = pct_threshold #Amount of variability desired to be explained by the principal components.  Set to match the results in book chapter and SAS code.  User can adjust this to a higher (>= 0.8) number but < 1.0
    #Load data
    theDataMatrix = data
    dataRowN = nrow(theDataMatrix)
    dataColN = ncol(theDataMatrix)
    #Center the data (center rows)
    theDataMatrixCentered = matrix(data = 0, nrow = dataRowN, ncol = dataColN)
    theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
    theDataMatrixCentered = t(theDataMatrixCentered_transposed)
    #trait data
    exp_design = trait
    rownames(exp_design) = rownames(trait)
    expDesignRowN = nrow(exp_design)
    expDesignColN = ncol(exp_design)
    myColNames = names(exp_design)
    #Compute correlation matrix
    theDataCor = cor(theDataMatrixCentered)
    #Obtain eigenvalues
    eigenData = eigen(theDataCor)
    eigenValues = eigenData$values
    ev_n = length(eigenValues)
    eigenVectorsMatrix = eigenData$vectors
    eigenValuesSum = sum(eigenValues)
    percents_PCs = eigenValues /eigenValuesSum
    #Merge experimental file and eigenvectors for n components
    my_counter_2 = 0
    my_sum_2 = 1
    for (i in ev_n:1){
        my_sum_2  = my_sum_2 - percents_PCs[i]
        if ((my_sum_2) <= pct_threshold ){
            my_counter_2 = my_counter_2 + 1
        }
    }
    if (my_counter_2 < 3){
        pc_n  = 3
    }else {
        pc_n = my_counter_2
    } #pc_n is the number of principal components to model
    pc_data_matrix = matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
    mycounter = 0
    for (i in 1:pc_n){
        for (j in 1:expDesignRowN){
            mycounter = mycounter + 1
            pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
        }
    }
    AAA = exp_design[rep(1:expDesignRowN,pc_n),]
    Data = cbind(AAA,pc_data_matrix)
    #Edit these variables according to your factors
    # Data$sex = as.factor(Data$sex)
    variables = c(colnames(exp_design))
    for (i in 1:length(variables)) {
        Data$variables[i] = as.factor(Data$variables[i])
    }
    #Mixed linear model
    op = options(warn = (-1))
    effects_n = expDesignColN + 1
    randomEffectsMatrix = matrix(data = 0, nrow = pc_n, ncol = effects_n)
    model.func = c()
    index = 1
    for (i in 1:length(variables)) {
        mod = paste("(1|", variables[i], ")", sep = "")
        model.func[index] = mod
        index = index + 1
    }
    function.mods = paste(model.func, collapse = " + ")
    for (i in 1:pc_n) {
        y = (((i - 1) * expDesignRowN) + 1)
        funct = paste("pc_data_matrix", function.mods, sep = " ~ ")
        Rm1ML = lme4::lmer(funct, Data[y:(((i - 1) * expDesignRowN) + expDesignRowN), ], REML = TRUE,
            control=lme4::lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"),verbose = FALSE, na.action = na.omit)
        randomEffects = Rm1ML
        randomEffectsMatrix[i, ] = c(unlist(lme4::VarCorr(Rm1ML)), resid = sigma(Rm1ML)^2)
    }
    effectsNames = c(names(lme4::getME(Rm1ML, "cnms")), "resid")
    #Standardize Variance
    randomEffectsMatrixStdze = matrix(data = 0, nrow = pc_n, ncol = effects_n)
    for (i in 1:pc_n){
        mySum = sum(randomEffectsMatrix[i,])
        for (j in 1:effects_n){
            randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum
        }
    }
    #Compute Weighted Proportions
    randomEffectsMatrixWtProp = matrix(data = 0, nrow = pc_n, ncol = effects_n)
    for (i in 1:pc_n){
        weight = eigenValues[i]/eigenValuesSum
        for (j in 1:effects_n){
            randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
        }
    }
    #Compute Weighted Ave Proportions
    randomEffectsSums = matrix(data = 0, nrow = 1, ncol = effects_n)
    randomEffectsSums = colSums(randomEffectsMatrixWtProp)
    totalSum = sum(randomEffectsSums)
    randomEffectsMatrixWtAveProp = matrix(data = 0, nrow = 1, ncol = effects_n)
    for (j in 1:effects_n){
        randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
    }
    bp = barplot(randomEffectsMatrixWtAveProp,  main = main, xlab = xlab, ylab = "Weighted average proportion variance", ylim= c(0,1.1),col = c("blue"), las=2)
    #replace the below code of "axis(1, at = bp, labels = effectsNames, xlab = "Effects", cex.axis = 0.5, las=2)" if you want rotate the x axis labels.
    axis(1, at = bp, labels = effectsNames, xlab = xlab, cex.axis = 0.8, las=2)
    #text(bp, par("usr")[3]-0.02, srt=45, adj=1, labels=effectsNames, xpd=TRUE, cex=0.8)
    values = randomEffectsMatrixWtAveProp
    new_values = round(values, 3)
    text(bp,randomEffectsMatrixWtAveProp,labels = new_values, pos=3, cex = 0.8) # place numbers on top of bars
}

#' Principal Variance Component Analysis 2
#'
#' This function perform PVCA with interaction.
#'
#' @param data Data matrix
#' @param trait Sample information used in the liner model
#' @param pct_threshold The threshold for PC Defaults to .5876
#' @param main The title in plot Defaults to "PVCA plot"
#' @param xlab The xlab in plot Defaults to "Effects"
#' @keywords PVCA
#' @return PVCA plot
#' @export
#' @examples
#' pvca2(data, trait)
#' @export

pvca2 <- function(data,trait,pct_threshold = .5876, main="PVCA plot", xlab="Effects"){
    pct_threshold = pct_threshold #Amount of variability desired to be explained by the principal components.  Set to match the results in book chapter and SAS code.  User can adjust this to a higher (>= 0.8) number but < 1.0
    #Load data
    theDataMatrix = data
    dataRowN = nrow(theDataMatrix)
    dataColN = ncol(theDataMatrix)
    #Center the data (center rows)
    theDataMatrixCentered = matrix(data = 0, nrow = dataRowN, ncol = dataColN)
    theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
    theDataMatrixCentered = t(theDataMatrixCentered_transposed)
    #Trait data
    exp_design = trait
    rownames(exp_design) = rownames(trait)
    expDesignRowN = nrow(exp_design)
    expDesignColN = ncol(exp_design)
    myColNames = names(exp_design)
    #Compute correlation matrix
    theDataCor = cor(theDataMatrixCentered)
    #Obtain eigenvalues
    eigenData = eigen(theDataCor)
    eigenValues = eigenData$values
    ev_n = length(eigenValues)
    eigenVectorsMatrix = eigenData$vectors
    eigenValuesSum = sum(eigenValues)
    percents_PCs = eigenValues /eigenValuesSum
    #Merge experimental file and eigenvectors for n components
    my_counter_2 = 0
    my_sum_2 = 1
    for (i in ev_n:1){
        my_sum_2  = my_sum_2 - percents_PCs[i]
        if ((my_sum_2) <= pct_threshold ){
            my_counter_2 = my_counter_2 + 1
        }
    }
    if (my_counter_2 < 3){
        pc_n  = 3
    }else {
        pc_n = my_counter_2
    } #pc_n is the number of principal components to model
    pc_data_matrix = matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
    mycounter = 0
    for (i in 1:pc_n){
        for (j in 1:expDesignRowN){
            mycounter = mycounter + 1
            pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
        }
    }
    AAA = exp_design[rep(1:expDesignRowN,pc_n),]
    Data = cbind(AAA,pc_data_matrix)
    #Edit these variables according to your factors
    # Data$sex = as.factor(Data$sex)
    variables = c(colnames(exp_design))
    for (i in 1:length(variables)) {
        Data$variables[i] = as.factor(Data$variables[i])
    }
    #Mixed linear model
    op = options(warn = (-1))
    effects_n = expDesignColN + choose(expDesignColN, 2) + 1
    randomEffectsMatrix = matrix(data = 0, nrow = pc_n, ncol = effects_n)
    model.func = c()
    index = 1
    for (i in 1:length(variables)) {
        mod = paste("(1|", variables[i], ")", sep = "")
        model.func[index] = mod
        index = index + 1
    }
    for (i in 1:(length(variables) - 1)) {
        for (j in (i + 1):length(variables)) {
            mod = paste("(1|", variables[i], ":", variables[j], ")", sep = "")
            model.func[index] = mod
            index = index + 1
        }
    }
    function.mods = paste(model.func, collapse = " + ")
    for (i in 1:pc_n) {
        y = (((i - 1) * expDesignRowN) + 1)
        funct = paste("pc_data_matrix", function.mods, sep = " ~ ")
        Rm1ML = lme4::lmer(funct, Data[y:(((i - 1) * expDesignRowN) + expDesignRowN), ], REML = TRUE,
            control=lme4::lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"),verbose = FALSE, na.action = na.omit)
        randomEffects = Rm1ML
        randomEffectsMatrix[i, ] = c(unlist(lme4::VarCorr(Rm1ML)), resid = sigma(Rm1ML)^2)
    }
    effectsNames = c(names(lme4::getME(Rm1ML, "cnms")), "resid")
    #Standardize Variance
    randomEffectsMatrixStdze = matrix(data = 0, nrow = pc_n, ncol = effects_n)
    for (i in 1:pc_n){
        mySum = sum(randomEffectsMatrix[i,])
        for (j in 1:effects_n){
            randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum
        }
    }
    #Compute Weighted Proportions
    randomEffectsMatrixWtProp = matrix(data = 0, nrow = pc_n, ncol = effects_n)
    for (i in 1:pc_n){
        weight = eigenValues[i]/eigenValuesSum
        for (j in 1:effects_n){
            randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
        }
    }
    #Compute Weighted Ave Proportions
    randomEffectsSums = matrix(data = 0, nrow = 1, ncol = effects_n)
    randomEffectsSums = colSums(randomEffectsMatrixWtProp)
    totalSum = sum(randomEffectsSums)
    randomEffectsMatrixWtAveProp = matrix(data = 0, nrow = 1, ncol = effects_n)
    for (j in 1:effects_n){
        randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
    }
    bp = barplot(randomEffectsMatrixWtAveProp,  main = main, xlab = xlab, ylab = "Weighted average proportion variance", ylim= c(0,1.1),col = c("blue"), las=2)
    #replace the below code of "axis(1, at = bp, labels = effectsNames, xlab = "Effects", cex.axis = 0.5, las=2)" if you want rotate the x axis labels.
    axis(1, at = bp, labels = effectsNames, xlab = xlab, cex.axis = 0.8, las=2)
    # text(bp, par("usr")[3]-0.02, srt = 45, adj = 1,labels = effectsNames, xpd = TRUE,cex=0.8)
    values = randomEffectsMatrixWtAveProp
    new_values = round(values , 3)
    text(bp,randomEffectsMatrixWtAveProp,labels = new_values, pos=3, cex = 0.8) #place numbers on top of bars
}

