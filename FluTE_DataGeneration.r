#############################
### FluTE Data Generation ###
#############################

# This script generates a training set and a test set for experimentation.
# For the details how to set up the FluTE model to your computer, the reader is referred to https://github.com/dlchao/FluTE
# You may also parallelize the FluTE model runs by using the parallel package in R (https://stat.ethz.ch/R-manual/R-patched/RHOME/library/parallel/html/parallel-package.html)
# You may also use the parallel version of the FluTE called mpiflute.

# Required R packages
# If they are not installed prior to loading, you need to install them first by using install.packages() function
# e.g. install.packages("lhs")
library(lhs)			# for maximin Latin hypercube sampling

################################
### Training Data Generation ###
################################

# Generate an empty matrix to store input parameter combinations
# Assume there are 30 input parameter combinations in the training set
tr = matrix(NA, nrow = 30, ncol = 7)	

# Generate input parameter combinations with 30 instances using maximin LHS
set.seed(1)
design = maximinLHS(n = 30, k = 7, dup = 5)
tr[, 1] = qunif(design[, 1], 1.2, 2.4) # R0
tr[, 2] = round(qunif(design[, 2], 0, 10000)) # antiviraldosesdaily
tr[, 3] = round(qunif(design[, 3], 0, 10000)) # vaccinedosesdaily
tr[, 4] = qunif(design[, 4], 0, 1) # isolation
tr[, 5] = qunif(design[, 5], 0, 1) # quarantine
tr[, 6] = qunif(design[, 6], 0, 1) # liberalleave
tr[, 7] = round(qunif(design[, 7], 0, 90)) # schoolclosuredays
tr = as.data.frame(tr)

colnames(tr) = c("R0", "antiviraldosesdaily", "vaccinedosesdaily", "isolation", "quarantine", "liberalleave", "schoolclosuredays")

# Set the other required parameters to run the FluTE
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

# Specify the number of replications for each parameter combination
nofrep = 20

results = NULL

# Run the model for each input parameter combination
for(i in 1:nrow(tr)){

	replresults = NULL

	# First, generate the required input files for the FluTE	
	for(j in 1:nofrep){
	
		sink("C:/path/to/save/input/file/config.txt")
		cat("label example")
		cat("\n")
		cat("datafile seattle")
		cat("\n")
		cat("seed ")
		cat(j)
		cat("\n")
		cat("R0 ")
		cat(as.numeric(tr[i, 1]))
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
		cat(as.numeric(tr[i, 2]))
		cat("\n")
		cat("vaccinedosesdaily ")
		cat(as.numeric(tr[i, 3]))
		cat("\n")
		cat("isolation ")
		cat(as.numeric(tr[i, 4]))
		cat("\n")
		cat("quarantine ")
		cat(as.numeric(tr[i, 5]))
		cat("\n")
		cat("liberalleave ")
		cat(as.numeric(tr[i, 6]))
		cat("\n")
		cat("schoolclosurepolicy ")
		cat(schoolclosurepolicy)
		cat("\n")
		cat("schoolclosuredays ")
		cat(as.numeric(tr[i, 7])) 
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
	
	# Store the replicated results with corresponding input parameter combinations
	results = rbind(results, cbind(tr[(rep(i, nofrep)), ], replresults))
	# Print the results
	print(results)
	
}

# Take the average of replications
tr = sapply(results, function(x) colMeans(matrix(x, nrow = nofrep)))
colnames(tr)[8] = "ar"
# Save the training set
write.csv(tr, "C:/path/to/save/training/data/flute_tr.csv") 

############################
### Test Data Generation ###
############################

# Generate an empty matrix to store input parameter combinations
# Assume there are 1000 input parameter combinations in the test set
test = matrix(NA, nrow = 1000, ncol = 7)	

# Generate input parameter combinations with 1000 instances using maximin LHS
set.seed(2)
design = maximinLHS(n = 1000, k = 7, dup = 5)
test[, 1] = qunif(design[, 1], 1.2, 2.4) # R0
test[, 2] = round(qunif(design[, 2], 0, 10000)) # antiviraldosesdaily
test[, 3] = round(qunif(design[, 3], 0, 10000)) # vaccinedosesdaily
test[, 4] = qunif(design[, 4], 0, 1) # isolation
test[, 5] = qunif(design[, 5], 0, 1) # quarantine
test[, 6] = qunif(design[, 6], 0, 1) # liberalleave
test[, 7] = round(qunif(design[, 7], 0, 90)) # schoolclosuredays
test = as.data.frame(test)

colnames(test) = c("R0", "antiviraldosesdaily", "vaccinedosesdaily", "isolation", "quarantine", "liberalleave", "schoolclosuredays")

# Set the other required parameters to run the FluTE
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

# Specify the number of replications for each parameter combination
nofrep = 20

results = NULL

# Run the model for each input parameter combination
for(i in 1:nrow(test)){

	replresults = NULL

	# First, generate the required input files for the FluTE	
	for(j in 1:nofrep){
	
		sink("C:/path/to/save/input/file/config.txt")
		cat("label example")
		cat("\n")
		cat("datafile seattle")
		cat("\n")
		cat("seed ")
		cat(j)
		cat("\n")
		cat("R0 ")
		cat(as.numeric(test[i, 1]))
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
		cat(as.numeric(test[i, 2]))
		cat("\n")
		cat("vaccinedosesdaily ")
		cat(as.numeric(test[i, 3]))
		cat("\n")
		cat("isolation ")
		cat(as.numeric(test[i, 4]))
		cat("\n")
		cat("quarantine ")
		cat(as.numeric(test[i, 5]))
		cat("\n")
		cat("liberalleave ")
		cat(as.numeric(test[i, 6]))
		cat("\n")
		cat("schoolclosurepolicy ")
		cat(schoolclosurepolicy)
		cat("\n")
		cat("schoolclosuredays ")
		cat(as.numeric(test[i, 7])) 
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
	
	# Store the replicated results with corresponding input parameter combinations
	results = rbind(results, cbind(test[(rep(i, nofrep)), ], replresults))
	# Print the results
	print(results)
	
}

# Take the average of replications
test = sapply(results, function(x) colMeans(matrix(x, nrow = nofrep)))
colnames(test)[8] = "ar"
# Save the test set
write.csv(test, "C:/path/to/save/test/data/flute_test.csv") 