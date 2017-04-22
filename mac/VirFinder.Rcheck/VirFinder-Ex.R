pkgname <- "VirFinder"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('VirFinder')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("VF.pred")
### * VF.pred

flush(stderr()); flush(stdout())

### Name: VF.pred
### Title: Identify viral sequences in a fasta file
### Aliases: VF.pred
### Keywords: default

### ** Examples


## (1) set the input fasta file name. 
library(VirFinder)
inFaFile <- system.file("data", "contigs.fa", package="VirFinder")

## (2) prediction
predResult <- VF.pred(inFaFile)
predResult

## sort sequences by p-value in ascending order
predResult[order(predResult$pvalue),]

## (3) predict for crAssphage
inFaFile <- system.file("data", "crAssphage.fa", package="VirFinder")
VF.pred(inFaFile)





cleanEx()
nameEx("VF.pred.user")
### * VF.pred.user

flush(stderr()); flush(stdout())

### Name: VF.pred.user
### Title: Identify viral sequences in a fasta file using user's trained
###   prediction model
### Aliases: VF.pred.user
### Keywords: customization

### ** Examples


## (1) specifiy the fasta files of the training contigs
#### (1.1) one for virus and one for prokaryotic hosts
trainFaFileHost <- system.file("data", "tara_host.fa", package="VirFinder")
trainFaFileVirus <- system.file("data", "tara_virus.fa", package="VirFinder")

#### (1.2) specify the directory where the trained model will be saved to, and the name of the model
userModDir <- file.path(find.package("VirFinder"))
userModName <- "modTara"

## (2) train the model using user's database
w <- 4  # the length of the k-tuple word
VF.trainModUser <- VF.train.user(trainFaFileHost, trainFaFileVirus, userModDir, 
userModName, w, equalSize=TRUE)

## (3) predict the contigs using the customized model
#### (3.1) specify the fasta file containing contigs for prediction
inFaFile <- system.file("data", "contigs.fa", package="VirFinder")

#### (3.2) prediction
predResultUser <- VF.pred.user(inFaFile, VF.trainModUser)
predResultUser

#### (3.3) sort sequences by p-value in ascending order
predResultUser[order(predResultUser$pvalue),]




cleanEx()
nameEx("VF.qvalue")
### * VF.qvalue

flush(stderr()); flush(stdout())

### Name: VF.qvalue
### Title: Estimate the false discovery rates (q-values) using p-values
### Aliases: VF.qvalue
### Keywords: FDR

### ** Examples


## (1) set the input fasta file name. 
library(VirFinder)
inFaFile <- system.file("data", "contigs.fa", package="VirFinder")

## (2) prediction
predResult <- VF.pred(inFaFile)
predResult

#### (2.1) sort sequences by p-value in ascending order
predResult[order(predResult$pvalue),]

#### (2.2) estimate q-values (false discovery rates) based on p-values
predResult$qvalue <- VF.qvalue(predResult$pvalue)

#### (2.3) sort sequences by q-value in ascending order
predResult[order(predResult$qvalue),]




cleanEx()
nameEx("VF.train.user")
### * VF.train.user

flush(stderr()); flush(stdout())

### Name: VF.train.user
### Title: Train virus prediction model using user's database
### Aliases: VF.train.user
### Keywords: customization

### ** Examples


## (1) train the model using user's database
#### (1.1) specifiy the fasta files of the training contigs, one for virus and one for prokaryotic hosts
trainFaFileHost <- system.file("data", "tara_host.fa", package="VirFinder")
trainFaFileVirus <- system.file("data", "tara_virus.fa", package="VirFinder")

#### (1.2) specify the directory where the trained model will be saved to, and the name of the model
userModDir <- file.path(find.package("VirFinder"))
userModName <- "modTara"

## (2) train the model using user's database
w <- 4  # the length of the k-tuple word
VF.trainModUser <- VF.train.user(trainFaFileHost, trainFaFileVirus, userModDir, 
userModName, w, equalSize=TRUE)

## (3) load the trained model based on user's database
modFile <- list.files(userModDir, userModName, full.names=TRUE)
load(modFile)




cleanEx()
nameEx("VF.trainMod8mer")
### * VF.trainMod8mer

flush(stderr()); flush(stdout())

### Name: VF.trainMod8mer
### Title: The prediction models trained using sequences of various lengths
###   from virus and host complete genomes.
### Aliases: VF.trainMod8mer
### Keywords: datasets

### ** Examples

data(VF.trainMod8mer)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
