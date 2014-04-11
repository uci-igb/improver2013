challenge <- 'subchallenge2' #'subchallenge3' # 'subchallenge1'

data.dir <- sprintf('../data/SBV_STC_%s/', challenge)
out.dir <- sprintf('../out/SBV_STC_%s/', challenge)

orthologs.file <- 'orthologs_used.txt'

rat.expression.file <- 'GEx_rat_train.txt'
rat.expression.file.test <- 'GEx_rat_test.txt'
human.expression.file <- 'GEx_human_train.txt'

rat.phospho.files <- c('Phospho_5_rat_train.txt', 'Phospho_25_rat_train.txt')
rat.phospho.files.test <- c('Phospho_5_rat_test.txt', 'Phospho_25_rat_test.txt')
human.phospho.files <- c('Phospho_5_human_train.txt', 'Phospho_25_human_train.txt')

rat.gene.sets <- c('gene_sets_fdr_rat_train.txt', 'gene_sets_nes_rat_train.txt')
rat.gene.sets.test <- c('gene_sets_fdr_rat_test.txt', 'gene_sets_nes_rat_test.txt')
human.gene.sets <- c('gene_sets_fdr_human_train.txt', 'gene_sets_nes_human_train.txt')

total.batches <- 8

root.dir <- dirname(getwd())
predictions.dir <- sprintf('%s/out/SBV_STC_%s/transformed/', root.dir, challenge)
test.file <- sprintf('%s/GEx_human_test.txt.predicted', predictions.dir)
train.file <- sprintf('%s/GEx_human_train.txt.predicted', predictions.dir)

test.samples.file <- sprintf('%s/GEx_rat_test.txt.sample_groups', out.dir)
train.samples.file <- sprintf('%s/GEx_rat_train.txt.sample_groups', out.dir)
