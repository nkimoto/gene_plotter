library(karyoploteR)

# args <- commandArgs(trailingOnly=TRUE)

selected_chr <- args[1]
reference <- args[2] # "Length_Silkbase.Sheet1.tsv"
annotation <- args[3] # "tblastx_HitTable_THlisted_trimABC_orthlist.txt"

pos.sb <- 0.1
pos.trinity <- 0.4
pos.bitscore <- 0.87
pos.evalue <- 0.95
size.chr <- 1.5
size.tick <- 0.8
size.sb <- 0.7
size.trinity <- 0.7
size.bitscore <- 0.7
size.evalue <- 0.7
size.left.label <- 0.6
color.sb <- "red"
color.trinity <- "black"
color.bitscore <- "blue"
color.evalue <- "purple"
fig.height.pdf <- 4.6
fig.width.pdf <- 20.0
fig.height.png <- 370
fig.width.png <- 1500

## selected_chr <- c("ChrZ")
## Write Reference genome
ref <- read.table(reference, sep="\t", row.names=1)
chr <- paste0("Chr", rownames(ref))
ref <- as.character(ref[,1])
start <- rep(1, length(chr))
end_tmp <- c()
end <- c()
for (i in strsplit(ref, "\\.\\.")){end_tmp <- append(end_tmp, i[2])}
for (i in strsplit(end_tmp, " \\(")){end <- append(end, i[1])}
end <- as.integer(end)
print(chr)
custom.genome <- toGRanges(data.frame(chr, start, end))

## Write Genes
genes <- read.table(annotation, header=TRUE, sep="\t", row.names=1)
genes["SB_chr"] <- sapply(strsplit(as.character(genes[["SB_chr"]]), "_"), function(x){x[2]})
genes[["SB_chr"]] <- sub("Chr1$", "ChrZ", genes[["SB_chr"]])

Trinity_ID <- rownames(genes)
SB_start <- genes[["SB_start"]]
SB_end <- genes[["SB_end"]]
print(genes[["SB_chr"]])
gene.symbols <- c(Trinity_ID)
granges <- makeGRangesFromDataFrame(genes,
                                    seqnames.field="SB_chr",
                                    start.field="SB_start",
                                    end.field="SB_end")
values(granges) <- genes[setdiff(colnames(genes), c("SB_chr", "SB_start", "SB_end"))]



pdf(paste0(selected_chr, ".pdf"), width=fig.width.pdf, height=fig.height.pdf)
kp <- plotKaryotype(genome=custom.genome, plot.type=4, cex=size.chr, chromosomes=selected_chr)
kpAddBaseNumbers(kp, tick.dist=1000000, tick.len=10, tick.col="red", cex=size.tick,
                 minor.tick.dist=100000, minor.tick.len=5, minor.tick.col="gray")

kpAddLabels(kp, "MetaData",
            label.margin=0.1,
            srt=90, pos=3, cex=0.8,
            data.panel=1)
kpAddLabels(kp, "SB_ID",
            r0=pos.sb + 0.15,
            r1=pos.sb + 0.15,
            cex=0.6,
            data.panel=1)
kpAddLabels(kp, "TRINITY_ID",
            r0=pos.trinity + 0.25,
            r1=pos.trinity + 0.25,
            cex=0.6,
            data.panel=1)
kpAddLabels(kp, "bitscore",
            r0=pos.bitscore + 0.025,
            r1=pos.bitscore + 0.025, cex=0.6,
            data.panel=1)
kpAddLabels(kp, "evalue",
            r0=pos.evalue + 0.025,
            r1=pos.evalue + 0.025,
            cex=0.6,
            data.panel=1)


kpPlotMarkers(kp,
              plot.type=2,
              data=granges,
              cex=0.7,
              labels=granges$SB_ID,
              r0=0,
              r1=pos.sb,
              label.color="red",
              marker.parts=c(0.2, 0.7, 0.1))
kpPlotMarkers(kp,
              plot.type=2,
              data=granges,
              cex=0.7, 
              labels=names(granges),
              r0=pos.trinity,
              r1=pos.trinity,
              label.color="black",
              marker.parts=c(0, 0, 0))
kpPlotMarkers(kp,
              plot.type=2,
              data=granges,
              cex=0.7, 
              labels=as.factor(granges$bitscore),
              r0=pos.bitscore,
              r1=pos.bitscore,
              marker.parts=c(0, 0, 0),
              label.color="blue")
kpPlotMarkers(kp,
              plot.type=2,
              data=granges,
              cex=0.7,
              labels=as.factor(granges$evalue),
              r0=pos.evalue,
              r1=pos.evalue,
              marker.parts=c(0, 0, 0),
              label.color="purple")
dev.off()
png(paste0(selected_chr, ".png"), width=fig.width.png, height=fig.height.png)
kp <- plotKaryotype(genome=custom.genome, plot.type=4, cex=size.chr, chromosomes=selected_chr)
kpAddBaseNumbers(kp, tick.dist=1000000, tick.len=10, tick.col="red", cex=size.tick,
                 minor.tick.dist=100000, minor.tick.len=5, minor.tick.col="gray")

kpAddLabels(kp, "MetaData",
            label.margin=0.1,
            srt=90, pos=3, cex=0.8,
            data.panel=1)
kpAddLabels(kp, "SB_ID",
            r0=pos.sb + 0.15,
            r1=pos.sb + 0.15,
            cex=0.6,
            data.panel=1)
kpAddLabels(kp, "TRINITY_ID",
            r0=pos.trinity + 0.25,
            r1=pos.trinity + 0.25,
            cex=0.6,
            data.panel=1)
kpAddLabels(kp, "bitscore",
            r0=pos.bitscore + 0.025,
            r1=pos.bitscore + 0.025, cex=0.6,
            data.panel=1)
kpAddLabels(kp, "evalue",
            r0=pos.evalue + 0.025,
            r1=pos.evalue + 0.025,
            cex=0.6,
            data.panel=1)


kpPlotMarkers(kp,
              plot.type=2,
              data=granges,
              cex=0.7,
              labels=granges$SB_ID,
              r0=0,
              r1=pos.sb,
              label.color="red",
              marker.parts=c(0.2, 0.7, 0.1))
kpPlotMarkers(kp,
              plot.type=2,
              data=granges,
              cex=0.7, 
              labels=names(granges),
              r0=pos.trinity,
              r1=pos.trinity,
              label.color="black",
              marker.parts=c(0, 0, 0))
kpPlotMarkers(kp,
              plot.type=2,
              data=granges,
              cex=0.7, 
              labels=as.factor(granges$bitscore),
              r0=pos.bitscore,
              r1=pos.bitscore,
              marker.parts=c(0, 0, 0),
              label.color="blue")
kpPlotMarkers(kp,
              plot.type=2,
              data=granges,
              cex=0.7,
              labels=as.factor(granges$evalue),
              r0=pos.evalue,
              r1=pos.evalue,
              marker.parts=c(0, 0, 0),
              label.color="purple")
dev.off()