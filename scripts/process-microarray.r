source('common.r')

# get phospho data first
combined.phospho.data <- function (phospho.files) {
    phospho.5.data <- read.delim(sprintf('%s/%s', data.dir, phospho.files[1]), header=T)
    phospho.25.data <- read.delim(sprintf('%s/%s', data.dir, phospho.files[2]), header=T)
    phospho.5.data[,1] <- paste(phospho.5.data[,1], "5", sep='_')
    phospho.25.data[,1] <- paste(phospho.25.data[,1], "25", sep='_')
    # for challenege2, the test data did not have the same columns so
    # we need to remove any (1) that aren't in both
    shared.colnames <- intersect(colnames(phospho.5.data), colnames(phospho.25.data))
    phospho.5.data <- phospho.5.data[,match(shared.colnames, colnames(phospho.5.data))]
    phospho.25.data <- phospho.25.data[,match(shared.colnames, colnames(phospho.25.data))]
    rbind(phospho.5.data, phospho.25.data)
}

rat.phospho.data <- combined.phospho.data(rat.phospho.files)
if (challenge != 'subchallenge1') {
    rat.phospho.data.test <- combined.phospho.data(rat.phospho.files.test)
}
human.phospho.data <- combined.phospho.data(human.phospho.files)

combined.genelist.data <- function (genelist.files) {
    genelist.fdr.data <- read.delim(sprintf('%s/%s', data.dir, genelist.files[1]), header=T)
    genelist.nes.data <- read.delim(sprintf('%s/%s', data.dir, genelist.files[2]), header=T)
    genelist.fdr.data[,1] <- paste(genelist.fdr.data[,1], "fdr", sep='_')
    genelist.nes.data[,1] <- paste(genelist.nes.data[,1], "nes", sep='_')
    # for challenge2, the test data did not have the same columns so
    # we need to remove any (1) that aren't in both
    shared.colnames <- intersect(colnames(genelist.fdr.data)[-1], colnames(genelist.nes.data)[-1])
    genelist.fdr.data <- genelist.fdr.data[,c(1, match(shared.colnames, colnames(genelist.fdr.data)))]
    genelist.nes.data <- genelist.nes.data[,c(1, match(shared.colnames, colnames(genelist.nes.data)))]
    colnames(genelist.nes.data) <- colnames(genelist.fdr.data)
    rbind(genelist.fdr.data, genelist.nes.data)
}

rat.genelist.data <- combined.genelist.data(rat.gene.sets)
rat.genelist.data.test <- combined.genelist.data(rat.gene.sets.test)
human.genelist.data <- combined.genelist.data(human.gene.sets)

# now process microarray and combine with phospho data
process.microarray.file <- function (expression.file, phospho.data=NULL, genelist.data=NULL, batch.offset=0) {
    expression.data <- read.delim(sprintf('%s/%s', data.dir, expression.file), header=T)
    feature.names <- as.character(expression.data[,1])
    batch.num <- gsub("^.*_batch([0-9])$", "\\1", colnames(expression.data)[-1])
    batch.features <- sapply(as.numeric(batch.num)+batch.offset, function(b) as.numeric(1:total.batches %in% b))
    output.data <- cbind(t(batch.features), t(expression.data[,-1]))
    feature.names <- c(paste("BatchNumber", 1:total.batches, sep=''), feature.names)
    write.table(feature.names, file=sprintf('%s/%s', out.dir, paste(expression.file, "features", sep='.')), row.names=F, quote=F, col.names=F)
    colnames(output.data) <- feature.names
    rownames(output.data) <- NULL
    dir.create(out.dir, showWarnings=F)

    clean.cols <- function (data) gsub("_rep[0-9]_batch[0-9]", "", colnames(data)[-1])
    clean.cols2 <- function (data) gsub("_batch[0-9]", "", colnames(data)[-1])
    # make sure the phospho and microarray data match
    if (!is.null(phospho.data) && is.null(genelist.data)) {
        output.phospho.data <- t(phospho.data[,-1])
        colnames(output.phospho.data) <- phospho.data[,1]
        phospho.cols <- clean.cols(phospho.data)
        expression.cols <- clean.cols(expression.data)

        # match on names alone, and combine all pairs of data into the same rows
        matches <- sapply(expression.cols, function (col.name) which(phospho.cols==col.name))
        expression.rows <- rep(1:length(expression.cols), times=sapply(matches, length))
        phospho.rows <- unlist(matches)

        output.data.full <- output.data[expression.rows,]
        output.data.cols <- expression.cols[expression.rows]
        output.phospho.data.full <- output.phospho.data[phospho.rows,]

        write.table(output.data.cols, file=sprintf('%s/%s', out.dir, paste(expression.file, "sample_groups", sep='.')), row.names=F, quote=F, col.names=F)
        write.table(output.data.full, file=sprintf('%s/%s', out.dir, expression.file), row.names=F, quote=F)
        write.table(output.phospho.data.full, file=sprintf('%s/%s', out.dir, paste(expression.file, "output", sep='.')), row.names=F, quote=F)
    }
    # make sure the genelist and microarray data match (and hence all match)
    # need to repeat genelist data for each sample replicate
    if (!is.null(genelist.data) && !is.null(phospho.data)) {
        output.phospho.data <- t(phospho.data[,-1])
        colnames(output.phospho.data) <- phospho.data[,1]
        phospho.cols <- clean.cols(phospho.data)
        expression.cols <- clean.cols(expression.data)

        # match on names alone, and combine all pairs of data into the same rows
        matches <- sapply(expression.cols, function (col.name) which(phospho.cols==col.name))
        expression.rows <- rep(1:length(expression.cols), times=sapply(matches, length))
        phospho.rows <- unlist(matches)

        output.data.full <- output.data[expression.rows,]
        output.data.cols <- expression.cols[expression.rows]
        output.phospho.data.full <- output.phospho.data[phospho.rows,]

        write.table(output.data.cols, file=sprintf('%s/%s', out.dir, paste(expression.file, "sample_groups", sep='.')), row.names=F, quote=F, col.names=F)
        write.table(output.data.full, file=sprintf('%s/%s', out.dir, expression.file), row.names=F, quote=F)
        write.table(output.phospho.data.full, file=sprintf('%s/%s', out.dir, paste(expression.file, "phospho", sep='.')), row.names=F, quote=F)

        output.genelist.data <- t(genelist.data[,-1])
        colnames(output.genelist.data) <- genelist.data[,1]
        genelist.cols <- clean.cols2(genelist.data)
        # match output.genelist to the stim type in output.data.cols

        print(output.data.cols)
        print(genelist.cols)
        matches <- match(output.data.cols, genelist.cols)

        print(matches)
        output.genelist.data.full <- output.genelist.data[matches,]
        output.genelist.data.nes <- output.genelist.data.full[,grep("_nes", colnames(output.genelist.data.full))]
        output.genelist.data.fdr <- output.genelist.data.full[,grep("_fdr", colnames(output.genelist.data.full))]

        write.table(output.genelist.data.full, file=sprintf('%s/%s', out.dir, paste(expression.file, "nes_and_fdr", sep='.')), row.names=F, quote=F)
        write.table(output.genelist.data.nes, file=sprintf('%s/%s', out.dir, paste(expression.file, "nes", sep='.')), row.names=F, quote=F)
        write.table(output.genelist.data.fdr, file=sprintf('%s/%s', out.dir, paste(expression.file, "fdr", sep='.')), row.names=F, quote=F)
    }
    if (is.null(phospho.data) && is.null(genelist.data)) {
        write.table(clean.cols(expression.data), file=sprintf('%s/%s', out.dir, paste(expression.file, "sample_groups", sep='.')), row.names=F, quote=F, col.names=F)
        write.table(output.data, file=sprintf('%s/%s', out.dir, expression.file), row.names=F, quote=F)
    }
}

if (challenge == 'subchallenge3') {
    process.microarray.file(rat.expression.file, rat.phospho.data, rat.genelist.data)
    process.microarray.file(human.expression.file, human.phospho.data, human.genelist.data, batch.offset=4)
    process.microarray.file(rat.expression.file.test, rat.phospho.data.test, rat.genelist.data.test)
} else {
    process.microarray.file(rat.expression.file, rat.phospho.data)
    process.microarray.file(human.expression.file, human.phospho.data, batch.offset=4)

    if (challenge == 'subchallenege2') {
        process.microarray.file(rat.expression.file.test, rat.genelist.data.test)
    } else {
        process.microarray.file(rat.expression.file.test)
    }
}
