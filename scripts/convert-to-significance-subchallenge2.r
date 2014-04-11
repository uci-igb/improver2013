##
# Read in the phosphorylation values
# 1) Training data converted to 0-1 by clipping between 0 to 4 and dividing by 4
# 2) Test data predicted by the model
#
# Predict 0-1 for changing, or not changing, across replicates as compared to DME
##
source("../opt/cybtpy/static/r-code/bayesreg.R")

library(gplots)
library(VGAM)

convert <- function (test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='', do.transform=NA, do.test='ttest', do.aggregate='max') {
    submission.out.dir <- sprintf('%s/%s', submission.out.dir, dir)
    dir.create(submission.out.dir)
    test.data <- read.table(test.file)
    train.data <- read.table(train.file)

    test.samples <- read.table(test.samples.file)$V1
    train.samples <- read.table(train.samples.file)$V1

    control.samples <- which(train.samples=="DME")
    control.data <- unique(train.data[control.samples,])

    output.data <- read.table(sprintf('%s/GEx_rat_train.txt.output', out.dir), header=T)
    output.cols <- colnames(output.data)

    colnames(control.data) <- output.cols
    colnames(test.data) <- output.cols

    raw.treatment.data <- list()
    for (samplename in levels(factor(test.samples))) {
        raw.treatment.data[[samplename]] <- t(rbind(control.data, test.data[which(test.samples==samplename),]))
    }

    fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }
    transform.r2z <- function(data) { apply(data, 1:2, function(r) fisher.r2z(min(r, 1-0.000000001))) }

    control.data <- transform.r2z(control.data)
    test.data <- transform.r2z(test.data)

    num.controls <- nrow(control.data)

    treatment.data <- list()
    for (samplename in levels(factor(test.samples))) {
        treatment.data[[samplename]] <- t(rbind(control.data, test.data[which(test.samples==samplename),]))
    }

    ttests <- sapply(treatment.data, function (df) {
        pval.transform <- function (p) log10(-log10(p)+1)
        df <- apply(df, 1:2, function (c) c+rnorm(1, sd=0.0001))
        
        df <- t(apply(df, 1, function (row) {
            if (is.na(do.transform)) {
                transform <- function (x) x
            } else {
                if (do.transform == 'logistic control')
                    transform <- function (x) plogis(x, location=max(control.data), scale=0.01)
                else if (do.transform == 'logistic 0.5')
                    transform <- function (x) plogis(x, location=0.5, scale=1)
            }
            transform(row)
        }))

        if (!is.na(do.test) & do.test == 'bayesian ttest') {
            num.controls <- nrow(control.data)
            p <- bayesT(as.data.frame(df), numC=num.controls, numE=ncol(df)-num.controls, winSize=11, conf=3, doMulttest=F)$pVal
        }

        if (is.na(do.test) | do.test == 'ttest') {
            ret <- apply(df, 1, function (row) {
                control.values <- row[1:num.controls]
                treatment.values <- row[-(1:num.controls)]
                if (is.na(do.test)) {
                    mean(treatment.values)
                } else {
                    pval.transform(t.test(treatment.values, control.values, var.equal=T, alternative="greater")$p.value)
                }
            })
        } else { # 'bayesian ttest'
            ret <- sapply(p, pval.transform)
        }
        ret
    }, simplify=T)

    ttests[which(is.infinite(ttests))] <- 1
    ttests[which(is.na(ttests))] <- 0

    rownames(ttests) <- colnames(output.data)

    png(sprintf('%s/%s-ttests.png', submission.out.dir, dir))
    heatmap.2(ttests, Colv=F, main='ttests', col=redgreen(100), trace='none')
    dev.off()

    results <- t(sapply(1:(nrow(ttests)/2), function (index) {
        pair.data <- ttests[c(index, index+16),]
        aggregate.func <- NA
        if (do.aggregate=='max')
            aggregate.func <- max
        else if (do.aggregate=='min')
            aggregate.func <- min
        else if (do.aggregate=='mean')
            aggregate.func <- mean
        
        apply(pair.data, 2, aggregate.func)
    }))
    rownames(results) <- gsub("_.*$", "", rownames(ttests)[1:(nrow(ttests)/2)])

    png(sprintf('%s/%s-results.png', submission.out.dir, dir))
    heatmap.2(results, Colv=F, main='results', col=redgreen(100), trace='none')
    dev.off()

    png(sprintf('%s/%s-results-hist.png', submission.out.dir, dir))
    hist.data <- hist(results, breaks=100)
    lines(rep(median(results), 2),c(0, max(hist.data$counts)), col='red', lwd=5)
    dev.off()
    
    submission.data <- (results-min(results))/max(results-min(results))
    rownames(submission.data) <- rownames(results)

    png(sprintf('%s/%s-submission.png', submission.out.dir, dir))
    hist.data <- hist(submission.data, breaks=10000, xlim=c(0,1))
    lines(rep(0.5, 2),c(0, max(hist.data$counts)), col='red', lwd=5)
    dev.off()

    png(sprintf('%s/%s-heatmap.png', submission.out.dir, dir))
    heatmap.2(submission.data, trace="none", Rowv=F, Colv=F, dendrogram="none", scale="none", col=greenred(100))
    dev.off()

    submission.data <- cbind(Phosphoprotein=rownames(submission.data), submission.data)
    write.table(submission.data, file=sprintf('%s/%s-submission.txt', submission.out.dir, dir), row.names=F, quote=F, sep='\t')
}

source('common.r')
submission.out.dir <- sprintf('%s/submissions', out.dir)

# this generates the submission.txt we submitted with:
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='original', do.transform=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='original-with-transform', do.transform='logistic control')
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='original-with-transform-0.5', do.transform='logistic 0.5')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='original-no-ttest-with-transform', do.transform='logistic control', do.test=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='original-no-ttest-with-transform-0.5', do.transform='logistic 0.5', do.test=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='original-no-ttest-with-transform-0.5-with-mean', do.transform='logistic 0.5', do.test=NA, do.aggregate='mean')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='original-no-ttest', do.transform=NA, do.test=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='original-with-bayes-ttest', do.transform=NA, do.test='bayesian ttest')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='original-with-min', do.transform=NA, do.aggregate="min")
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='original-with-mean', do.transform=NA, do.aggregate="mean")

source('common-use-rat.r')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='rat-only', do.transform=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='rat-only-with-transform', do.transform='logistic control')
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='rat-only-with-transform-0.5', do.transform='logistic 0.5')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='rat-only-no-ttest-with-transform', do.transform='logistic control', do.test=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='rat-only-no-ttest-with-transform-0.5', do.transform='logistic 0.5', do.test=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='rat-only-no-ttest-with-transform-0.5-with-mean', do.transform='logistic 0.5', do.test=NA, do.aggregate='mean')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='rat-only-no-ttest', do.transform=NA, do.test=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='rat-only-with-bayes-ttest', do.transform=NA, do.test='bayesian ttest')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='rat-only-with-min', do.transform=NA, do.aggregate="min")
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='rat-only-with-mean', do.transform=NA, do.aggregate="mean")

source('common-use-GEx.r')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='use-GEx', do.transform=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='use-GEx-with-transform', do.transform='logistic control')
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='use-GEx-with-transform-0.5', do.transform='logistic 0.5')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='use-GEx-no-ttest-with-transform', do.transform='logistic control', do.test=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='use-GEx-no-ttest-with-transform-0.5', do.transform='logistic 0.5', do.test=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='use-GEx-no-ttest-with-transform-0.5-with-mean', do.transform='logistic 0.5', do.test=NA, do.aggregate='mean')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='use-GEx-no-ttest', do.transform=NA, do.test=NA)
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='use-GEx-with-bayes-ttest', do.transform=NA, do.test='bayesian ttest')

convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='use-GEx-with-min', do.transform=NA, do.aggregate="min")
convert(test.file, test.samples.file, train.file, train.samples.file, submission.out.dir, dir='use-GEx-with-mean', do.transform=NA, do.aggregate="mean")


