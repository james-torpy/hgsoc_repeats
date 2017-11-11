library("rtracklayer")
library("GenomicRanges")

cGenes <- readRDS("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq//Robjects/exp2/c_RepeatGenes.rds")
crefDir <- paste0("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/refs/c")
system(paste0("mkdir -p ", crefDir))

custom3Genes <- readRDS("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/Robjects/exp2/custom3_RepeatGenes.rds")
custom3refDir <- paste0("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/refs/custom3")
system(paste0("mkdir -p ", custom3refDir))


cGTF <- lapply(cGenes, function(x) {
	nam <- gsub("/", "_", names(x))
	i=1
	gr <- lapply(x, function(y) {
		i=1
		strand(y) <- "*"
		reduce(y)
		y$ID <- seq(1, length(y))
		export(y, paste0("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/refs/c/", nam[i], ".gff"), "gff3")
		i <<- i+1
		return(y)
	})
})



custom3GTF <- lapply(custom3Genes, function(x) {
	i=1
	gr <- lapply(x, function(y) {
		nam <- gsub("/", "_", names(x)[i])
		strand(y) <- "*"
		y <- reduce(y)
		y$ID <- nam
		i <<- i+1
		return(y)
	})
	return(gr)
	j <<- j+1
})

custom3finalGTF <- do.call(getMethod(c, "GenomicRanges"),do.call("c", custom3GTF))
export(custom3finalGTF, paste0("/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/refs/custom3/custom3rep.gff"), "gff3")


# bash lines to change 'sequence_feature' to 'exon' in gffs:
#dir=`pwd`
#for f in custom3*; do echo $f; id=`echo $f | sed 's/.gff//'`; echo $dir/$id.final.gff; cat $f | sed 's/sequence_feature/exon/g' | sed 's/*/./g' > $dir/$id.final.gff; done

