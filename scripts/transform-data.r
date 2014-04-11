source('common.r')

for (file.base in c('rat_train', 'rat_test', 'human_train')) {
    output.file <- sprintf('%s/GEx_%s.txt.output', out.dir, file.base)
    transformed.file <- sprintf('%s/GEx_%s.txt.output', predictions.dir, file.base)
    data <- read.table(output.file, header=T)
    data <- apply(data, 1:2, function (c) min(max(0, c), 4)/4)
    write.table(data, file=transformed.file, col.names=F, row.names=F, quote=F, sep='\t')
}

