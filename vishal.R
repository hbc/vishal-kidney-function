library(readr)
fn = "data/vishal-7days.txt"
t = read_delim(fn, "\t")
symbols = unique(subset(t, padj < 0.1)$mgi_symbol)
write.table(unique(subset(t, padj < 0.1)$mgi_symbol), file="symbols.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
# convert symbols to U133 identifiers

library(biomaRt)
mouse = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
human = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
conversions = getLDS(attributes=c("ensembl_gene_id", "mgi_symbol"),
    attributesL=c("ensembl_gene_id", "hgnc_symbol", "affy_hg_u133a"),
    mart=mouse, martL=human)
autism = unique(subset(conversions, HGNC.symbol %in% autism$V1)$MGI.symbol)
retard = unique(subset(conversions, HGNC.symbol %in% retard$V1)$MGI.symbol)
seizure = unique(subset(conversions, HGNC.symbol %in% seizure$V1)$MGI.symbol)
#  affy_hg_u133a    is the array

up_symbols = unique(subset(t, padj < 0.1 & log2FoldChange > 0)$mgi_symbol)
up_u133 = subset(conversions, MGI.symbol %in% up_symbols)$Affy.HG.U133A.probeset
up_u133 = unique(up_u133)
down_symbols = unique(subset(t, padj < 0.1 & log2FoldChange < 0)$mgi_symbol)
down_u133 = subset(conversions, MGI.symbol %in% down_symbols)$Affy.HG.U133A.probeset
down_u133 = unique(down_u133)

write.table(up_u133, file="up_u133a.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(down_u133, file="down_u133a.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(up_symbols, file="up_symbols.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
write.table(down_symbols, file="down_symbols.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)

expressed = unique(subset(t, baseMean > 10)$mgi_symbol)
write.table(expressed, file="expressed_symbols.txt", quote=FALSE,
            row.names=FALSE, col.names=FALSE)
