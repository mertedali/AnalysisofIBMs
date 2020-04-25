######################
### FluTE Analysis ###
######################

# This script takes a training set and a test set as inputs.
# Then, the remaining part simply executes the Algorithm 1 presented in the manuscript.
# For the details how to set up FluTE model to your computer, the reader is referred to https://github.com/dlchao/FluTE
# You may also parallelize the FluTE model runs by using the parallel package in R (https://stat.ethz.ch/R-manual/R-patched/RHOME/library/parallel/html/parallel-package.html)
# You may also use the parallel version of the FluTE called mpiflute.

# Required R packages
# If they are not installed prior to loading, you need to install them first by using install.packages() function
# e.g. install.packages("randomForest")
library(randomForest)	# for random forest metamodeling	
library(e1071)			# for hyperparameter tuning of random forest
library(inTrees)		# for obtaining rule list from the random forest metamodel
library(gurobi)			# for solving the optimization problem given in Equation 4 in the manuscript
						# please note that gurobi is a commercial solver. however, it provides free one-year license for academic purposes

##########################################
### Active Learning with Random Forest ###
##########################################

# Define the error function (Root Mean Square Error in our case) to quantify metamodel accuracy
rmse = function(real, pred){
	error = pred - real
	return(sqrt(mean(error^2)))
}

# Read your training set
tr = read.csv("C:/path/to/training/data/flute_tr.csv", sep = ";")
# You may also use the readily available training data published on GitHub:
# tr = read.csv("https://raw.githubusercontent.com/mertedali/AnalysisofIBMs/master/flute_tr.csv", sep = ";")
# Read your test set
test = read.csv("C:/path/to/test/data/flute_test.csv", sep = ";")
# You may also use the readily available test data published on GitHub:
# test = read.csv("https://raw.githubusercontent.com/mertedali/AnalysisofIBMs/master/flute_test.csv", sep = ";")

# Warning: If you use the data given in GitHub, do not forget to take averages of runs.
# tr = sapply(tr, function(x) colMeans(matrix(x, nrow = nofrep)))
# test = sapply(test, function(x) colMeans(matrix(x, nrow = nofrep)))

# Make sure that both the training and test data are in the form of data frame
tr = as.data.frame(tr)
test = as.data.frame(test)

# Set seet for reproducibility
set.seed(1)
# Tune the hyperparameters of the random forest metamodel
tuning = tune.randomForest(ar ~ ., data = tr, ntree = c(100, 300, 500, 700), mtry = c(1, 3, 5, 7), tunecontrol = tune.control(sampling = "cross", cross = nrow(tr), error.fun = rmse))
# Train the random forest metamodel with the optimum hyperparameter combination
rf = randomForest(ar ~ ., data = tr, ntree = as.numeric(tuning$best.parameters$ntree), mtry = as.numeric(tuning$best.parameters$mtry))
# Generate predictions for test data instances
pred = predict(rf, newdata = test[, 1:7])
# Determine the accuracy of the metamodel on the test set
rmse(test$ar, pred)

# Generate a matrix to save accuracy on the test set throughout iterations
acc = matrix(NA, nrow = 15, ncol = 1)
acc[1, 1] = rmse(test$ar, pred)

# Generate an unlabeled dataset with 500 instances

unl = matrix(NA, nrow = 500, ncol = 7)

set.seed(1)
design = maximinLHS(n = 500, k = 7, dup = 5)
unl[, 1] = qunif(design[, 1], 1.2, 2.4) # R0
unl[, 2] = round(qunif(design[, 2], 0, 10000)) # antiviraldosesdaily
unl[, 3] = round(qunif(design[, 3], 0, 10000)) # vaccinedosesdaily
unl[, 4] = qunif(design[, 4], 0, 1) # isolation
unl[, 5] = qunif(design[, 5], 0, 1) # quarantine
unl[, 6] = qunif(design[, 6], 0, 1) # liberalleave
unl[, 7] = round(qunif(design[, 7], 0, 90)) # schoolclosuredays
unl = as.data.frame(unl)

colnames(unl) = c("R0", "antiviraldosesdaily", "vaccinedosesdaily", "isolation", "quarantine", "liberalleave", "schoolclosuredays")

# Set the other required parameters of the FluTE
# Note that these are kept constant

seedset = seq(1, 20, 1)
runlength = 180
seedinfected = 10
vaccinedata = "0 0.4 0.4 0.67 0 0 0 0 0 0 0"
vaccinepriorities = "1 1 1 1 1 1 1 1 1 1 1 1 1"
AVEs = 0.3
AVEi = 0.62
AVEp = 0.60
antiviralpolicy = "treatmentonly"
reactivestrategy = "mass"
schoolclosurepolicy = "all"
vaccinationfraction = 1
responsethreshold = 0 
responsedelay = -1

nofrep = 20

# Start the active learning loop (14 iterations)
for (iter in 1:14){

	# Predict the outputs of unlabeled dataset instances
	predunl = predict(rf, newdata = unl, predict.all = TRUE)$individual

	# Calculate the uncertainty of predicitions (Step 3 in Algorithm 1)
	dev = NULL
	for (ii in 1:nrow(predunl)){
		dev = rbind(dev, sd(predunl[ii, ]))
	}
	# Determine the indices of unlabeled instances with the top highest 5 uncertainty values
	idx = which(dev >= sort(dev, decreasing = TRUE)[5], arr.ind = TRUE)[, 1][1:5]
	
	# Select those instances (Step 4 in Algorithm 1)
	# Here, run is equivalent to U'
	run = unl[idx, ]
	run = cbind(run, NA)
	colnames(run)[8] = c("ar")
	
	# Remove the selected unlabeled samples from the pool (Step 5 in Algorithm 1)
	# U <- U \ U'
	unl = unl[-idx, ]
	
	results = NULL
	for(i in 1:nrow(run)){
	
		replresults = NULL

		for(j in 1:nofrep){
		
			sink("C:/path/to/save/input/file/config.txt")
			cat("label example-minimal")
			cat("\n")
			cat("datafile seattle")
			cat("\n")
			cat("seed ")
			cat(j)
			cat("\n")
			cat("R0 ")
			cat(as.numeric(run[i, 1]))
			cat("\n")
			cat("vaccinationfraction ")
			cat(vaccinationfraction)
			cat("\n")
			cat("responsethreshhold ")
			cat(responsethreshold)
			cat("\n")
			cat("responsedelay ")
			cat(responsedelay)
			cat("\n")
			cat("reactivestrategy ")
			cat(noquote("mass"))
			cat("\n")
			cat("antiviraldosesdaily ")
			cat(as.numeric(run[i, 2]))
			cat("\n")
			cat("vaccinedosesdaily ")
			cat(as.numeric(run[i, 3]))
			cat("\n")
			cat("isolation ")
			cat(as.numeric(run[i, 4]))
			cat("\n")
			cat("quarantine ")
			cat(as.numeric(run[i, 5]))
			cat("\n")
			cat("liberalleave ")
			cat(as.numeric(run[i, 6]))
			cat("\n")
			cat("schoolclosurepolicy ")
			cat(schoolclosurepolicy)
			cat("\n")
			cat("schoolclosuredays ")
			cat(as.numeric(run[i, 7])) 
			cat("\n")
			cat("antiviralpolicy ")
			cat(noquote(antiviralpolicy))
			cat("\n")
			cat("AVEs ")
			cat(AVEs)
			cat("\n")
			cat("AVEi ")
			cat(AVEi)
			cat("\n")
			cat("AVEp ")
			cat(AVEp)
			cat("\n")
			cat("runlength ")
			cat(as.numeric(runlength))
			cat("\n")
			cat("seedinfected ")
			cat(seedinfected)
			cat("\n")
			cat("vaccinedata ")
			cat(noquote(vaccinedata))
			cat("\n")
			cat("vaccinepriorities ")
			cat(vaccinepriorities)
			sink()
			
			# Set your working directory to the folder where the FluTE model is located at
			setwd("C:/path/to/FluTE/directory/FluTE-master/")
			# Run the model with config file which includes the required parameter settings
			system("cmd.exe", input = "flute config.txt")
			
			# Read the first 100 lines of the output file
			outputfile = readLines("C:/path/to/FluTE/directory/FluTE-master/Summary0", 100)
			# Find the line where the output that we are looking for is located at
			arindex = grep("Total symptomatic attack rate:", outputfile)
			# Store the replicated outputs
			replresults = rbind(replresults, as.numeric(unlist(strsplit(outputfile[arindex], ":"))[2]))

		}

		# Add the FluTE model outputs to U'
		run[i, "ar"] = mean(replresults) 
		
	}
	
	# Expand the training set (Step 7 in Algorithm 1)
	# L <- L U U'
	tr = rbind(tr, run)
	
	# Retrain the metamodel (Step 8 in Algorithm 1)
	set.seed(100 + j)
	tuning = tune.randomForest(ar ~ ., data = tr, ntree = c(100, 300, 500, 700), mtry = c(1, 3, 5, 7), tunecontrol = tune.control(sampling = "cross", cross = 10, nrepeat = 10, error.fun = rmse))
	rf = randomForest(ar ~ ., data = tr, ntree = as.numeric(tuning$best.parameters$ntree), mtry = as.numeric(tuning$best.parameters$mtry))
	# Predict the outputs of test instances with the newly trained metamodel
	pred = predict(rf, newdata = test[, 1:7], type = "response")
	# Calculate the accuracy and store it
	acc[iter + 1, 1] = rmse(test$ar, pred)
	# Print the matrix which stores accuracy values
	print(acc)
	
}

# Save the expanded training set
write.csv(tr, "C:/path/to/save/training/set/tr.csv")

#######################
### Rule Extraction ###
#######################

# Generate a rule list (i.e., RL in the manuscript)
ruleExec = extractRules(RF2List(rf), tr[, 1:7], ntree = rf$ntree, maxdepth = 100) 
ruleMetric = getRuleMetric(ruleExec, tr[, 1:7], tr[, 8])
readableRules = presentRules(ruleMetric, colnames(tr[, 1:7]))

# Generate matrix A (Algorithm 2)
finalbinarym = NULL
cost = NULL

for(i in 1:dim(readableRules)[1]){
	l = as.character(readableRules[i, "condition"])
	termval = as.numeric(readableRules[i, "pred"])
	l = unlist(strsplit(l, "&"))
	evl = NULL
	for(j in 1:length(l)){
		evl = cbind(evl, eval(parse(text = paste("tr$", l[j]))))
	}
	evlfinal = NULL
	for(k in 1:dim(evl)[1]){
		evlfinal = rbind(evlfinal, all(evl[k, ]))
	}
	finalbinarym = cbind(finalbinarym, as.numeric(evlfinal))
	cost = c(cost, rmse(tr[which(as.numeric(evlfinal) == 1), "ar"], rep(termval, length(which(as.numeric(evlfinal) == 1)))))
}

colsums = colSums(finalbinarym)

### Generate a set of unlabeled data to calculate coverage approximately

n = 200000
covdata = matrix(NA, nrow = n, ncol = 7)
	
set.seed(1)
covdata[, 1] = runif(n, 1.2, 2.4) # R0
covdata[, 2] = round(runif(n, 0, 10000)) # antiviraldosesdaily
covdata[, 3] = round(runif(n, 0, 10000)) # vaccinedosesdaily
covdata[, 4] = runif(n, 0, 1) # isolation
covdata[, 5] = runif(n, 0, 1) # quarantine
covdata[, 6] = runif(n, 0, 1) # liberalleave
covdata[, 7] = round(runif(n, 0, 90)) # schoolclosuredays
covdata = as.data.frame(covdata)

colnames(covdata) = c("R0", "antiviraldosesdaily", "vaccinedosesdaily", "isolation", "quarantine", "liberalleave", "schoolclosuredays")

lambda = seq(0, 1, 0.05)

results = NULL
for(lam in 2:length(lambda)){

	model = list()

	model$A          = finalbinarym
	model$obj		 = cost / max(cost) + lambda[lam] * (1 / (colsums / max(colsums)))
	model$modelsense = "min"
	model$rhs        = rep(1, dim(tr)[1])
	model$sense      = rep("=", dim(tr)[1])
	model$vtype      = rep("B", ncol(finalbinarym))

	result = gurobi(model)
	
	# We set "result$x >= 0.5" since for some cases Gurobi fails to obtain optimum integer solution
	extcondset = readableRules[which(result$x >= 0.5), ]
	
	coverageidx = NULL
	for(i in 1:dim(extcondset)[1]){
		l = as.character(extcondset[i, "condition"])
		l = unlist(strsplit(l, "&"))
		evl = NULL
		for(j in 1:length(l)){
			evl = cbind(evl, eval(parse(text = paste("covdata$", l[j]))))
		}
		evlfinal = NULL
		for(k in 1:dim(evl)[1]){
			evlfinal = rbind(evlfinal, all(evl[k, ]))
		}
		coverageidx = c(coverageidx, as.numeric(rownames(covdata[which(evlfinal == TRUE), ])))
	}
	
	results = rbind(results, c(lambda[lam], length(unique(coverageidx))/nrow(covdata), nrow(extcondset), sum(sqrt(as.numeric(extcondset[, "err"])))))
	
}

results = as.data.frame(results)
colnames(results) = c("lambda", "cov", "nr", "te")

# Function for score calculation to perform weighted rule selection (Equation 5)
score = function(d, w1, w2, w3){
	dtemp = d
	dtemp$cov = (dtemp$cov - min(dtemp$cov))/(max(dtemp$cov) - min(dtemp$cov))
	dtemp$nr = (dtemp$nr - min(dtemp$nr))/(max(dtemp$nr) - min(dtemp$nr))
	dtemp$te = (dtemp$te - min(dtemp$te))/(max(dtemp$te) - min(dtemp$te))
	idx = which.max(w1 * dtemp$cov - w2 * dtemp$nr - w3 * dtemp$te)
	lst = list("opt" = d[idx, ], "vals" = cbind(d, w1 * dtemp$cov - w2 * dtemp$nr - w3 * dtemp$te))
	return(lst)
}

# Select the lambda value which generates the highest score
optlambda = score(results, 0.33, 0.33, 0.33)$opt

# Solve the rule extraction problem with the selected lambda value 

lambda = optlambda

time1 = proc.time()

library(gurobi)

model = list()

model$A          = finalbinarym
model$obj		 = cost / max(cost) + lambda * (1 / (colsums / max(colsums)))
model$modelsense = "min"
model$rhs        = rep(1, dim(tr)[1])
model$sense      = rep("=", dim(tr)[1])
model$vtype      = rep("B", ncol(finalbinarym))

result = gurobi(model)

print(result$objval)
print(result$x)

finalbinarym[, which(result$x == 1)]
rowSums(finalbinarym[, which(result$x == 1)])
colSums(finalbinarym[, which(result$x == 1)])

extcondset = readableRules[which(result$x == 1), ]
sum(as.numeric(extcondset[, "freq"]))
sum(as.numeric(extcondset[, "err"]))

extcondset[, c("condition", "pred", "err")]
