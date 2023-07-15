# To run:
#library(bigWig)
#source('TSSinference.R')
#gene.exon1 = parse.bed.exon1('gencode.hg38.firstExon.bed')
#potential.tss2 = TSSinference(gene.exon1, 'H9_HDACi_merged_minus_PE2.bigWig', 'H9_HDACi_merged_plus_PE2.bigWig')

# Running line by line: 
# bed=gene.exon1
# bw.plus='H9_HDACi_merged_minus_PE2.bigWig'
# bw.minus='H9_HDACi_merged_plus_PE2.bigWig'
# tssWin = 100
# top.num.peaks = 20
# low.limit.tss.counts = 3
# denWin = 100

TSSinference <- function(bed, bw.plus, bw.minus, tssWin = 100, top.num.peaks = 20, low.limit.tss.factor = 3, denWin = 100, densityFilter = FALSE) {
# When looking for upstream antisense TSSs: flip bw.plus/bw.minus files

# make bigger window for finding TSS positions since some might not be annotated
    bed[,2] = bed[,2] - tssWin
    bed[,3] = bed[,3] + tssWin
# load bigwigs
    bwPlus = load.bigWig(bw.plus)
    bwMinus = load.bigWig(bw.minus)

# Change low.limit.tss.counts to scale on user-provided factor based on lowest read count in the bigWigs
    low.limit.tss.counts = low.limit.tss.factor*bwPlus$min

# Creates a vector for each gene/bed entry in which (takes a few seconds)
# each entry in 'bed' (i think) becomes a position in the [first exon?] and each position's value
# is either "0" or some number corresponding to the sum of the number of reads that start there
    x.beds = bed6.step.bpQuery.bigWig(bwPlus, bwMinus, bed, abs.value = TRUE, step = 1, with.attributes = TRUE)

# Applies top.num() [see below] to the previous output (takes a few seconds)
# Creates a list of positions for each bed entry of the "local peaks" of read pile-ups. 
    v = lapply(x.beds, top.num, top.num.peaks = top.num.peaks, low.limit.tss.counts = low.limit.tss.counts)
# Transpose object 'v' and convert to data frame so each item in 'v' becomes its own row 
# Number of columns is determined by "top.num.peaks" which defaults to 20.
# `w` becomes an index for which genes have peaks   
    w = as.data.frame(t(sapply(v, "[", i = seq_len(top.num.peaks))))

                                        #now make four data frames: all NA, top.num.peaks - 1 NAs, top.num.peaks -2 NAs, all others
                                        #this is the defaulting to the upstream most TSS instances

# Making a bed of all genes with NO peaks
    no.peaks = x.beds[is.na(w[,1])]
    no.peaks.bed = bed[is.na(w[,1]),]
# Gene start sites are returned to their original annotated values based on whether they're plus- or minus-stranded. 
    no.peaks.bed[,2][no.peaks.bed$strand == '+'] = no.peaks.bed[,2][no.peaks.bed$strand == '+'] + tssWin
    no.peaks.bed[,3][no.peaks.bed$strand == '-'] = no.peaks.bed[,3][no.peaks.bed$strand == '-'] - tssWin
# And then adjusted so the start/end reflects only 1 bp--the position of the beginning of the gene. 
    no.peaks.bed[,3][no.peaks.bed$strand == '+'] =  no.peaks.bed[,2][no.peaks.bed$strand == '+'] + 1
    no.peaks.bed[,2][no.peaks.bed$strand == '-'] =  no.peaks.bed[,3][no.peaks.bed$strand == '-'] - 1 
    
# Because this is the "no peaks" data, add columns but have them all be NA
    no.peaks.bed$height <- NA
    no.peaks.bed$up <- NA
    no.peaks.bed$down <- NA
#    print('no.peaks.bed')
#    print(head(no.peaks.bed))

# Making a bed of all genes with only 1 peak
    one.peak = w[!is.na(w[,1]) & is.na(w[,2]),]
    one.peak.bed = bed[!is.na(w[,1]) & is.na(w[,2]),]
    
    one.peak.bed.transform = one.peak.bed
# This is making a new vector with all the peak coordinates
# Change "end" position to position of sole peak    
    one.peak.bed.transform[,3] = one.peak.bed[,2] + one.peak[,1]
# Change "start" position to be peak-1, irrespective of strand
    one.peak.bed.transform[,2] = one.peak.bed.transform[,3] - 1
    
# Upstream.bed is stored in bed_utils.R in the bigWig repo, and this is making a bed with location of the upstream window.
    upstream.x = upstream.bed(one.peak.bed.transform, denWin)
# Makes a list of upstream density based on the upstream windows (upstream.x) for every gene.
    up.density = bed6.region.bpQuery.bigWig(bwPlus, bwMinus, upstream.x, op = "avg", abs.value = TRUE, gap.value = 0, bwMap = NULL)

# Makes bed of locations of downstream windows
    downstream.x = downstream.bed(one.peak.bed.transform, denWin)
# Uses downstream windows to calculate downstream densities for each gene.
    down.density = bed6.region.bpQuery.bigWig(bwPlus, bwMinus, downstream.x, op = "avg", abs.value = TRUE, gap.value = 0, bwMap = NULL)
# Uses the bigwigs to calculate the total number of reads at the peak for these one-peak genes
    peak.height= bed6.region.bpQuery.bigWig(bwPlus, bwMinus, one.peak.bed.transform, op = "sum", abs.value = TRUE, gap.value = 0, bwMap = NULL)  

    #!!!!! Run only if densityFilter = TRUE !!!! Also change for all density comparisons!!!
    if (densityFilter == TRUE) {
        # Change TSS start locations back to original positions (defaulting to upstream-most annotated position)
        one.peak.bed[,2][one.peak.bed$strand == '+'] = one.peak.bed[,2][one.peak.bed$strand == '+'] + tssWin
        one.peak.bed[,3][one.peak.bed$strand == '-'] = one.peak.bed[,3][one.peak.bed$strand == '-'] - tssWin

        # If down > up density, reassign TSS to match location of the one peak
        one.peak.bed[,2][down.density > up.density & one.peak.bed$strand == '+'] <- one.peak.bed.transform[,2][down.density > up.density & one.peak.bed$strand == '+']
        one.peak.bed[,3][down.density > up.density & one.peak.bed$strand == '-'] <- one.peak.bed.transform[,3][down.density > up.density & one.peak.bed$strand == '-']

        # Adjust start/ends so TSS has length 1bp
        one.peak.bed[,3][one.peak.bed$strand == '+'] = one.peak.bed[,2][one.peak.bed$strand == '+'] + 1
        one.peak.bed[,2][one.peak.bed$strand == '-'] = one.peak.bed[,3][one.peak.bed$strand == '-'] - 1

        # Add height, up density, and down density information to the table
        # Will have NAs for genes that didn't meet the (down > up) density criteria if this option is chosen
        one.peak.bed$height[down.density > up.density] = peak.height[down.density > up.density]
        one.peak.bed$up[down.density > up.density] = up.density[down.density > up.density]
        one.peak.bed$down[down.density > up.density] = down.density[down.density > up.density]
    } else {
        one.peak.bed = one.peak.bed.transform
        one.peak.bed$height = peak.height
        one.peak.bed$up = up.density
        one.peak.bed$down = down.density
    }

#    print('one.peak.bed')
#    print(head(one.peak.bed, 20))

    # two.peaks is actually all genes with 2+ peaks
    two.peaks = w[!is.na(w[,3]),]
    two.peaks.bed = bed[!is.na(w[,3]),] # gets annotated exon 1 boundaries for each of these genes
    # create table where each gene is followed by its peak location columns
    all.two.peaks.bed = cbind(two.peaks.bed, two.peaks)
    
    #this is complicated, but quick
    # Create a new column of all the relative locations of the top peaks separated by colon
    all.two.peaks.bed$newcol <- apply(all.two.peaks.bed[, c(7:(6+top.num.peaks))], 1,
                      function(i){ paste(na.omit(i), collapse = ":") })
    
    # This stores the "newcol" information with the individual locations separate; one row? per gene
    spikes <- strsplit(all.two.peaks.bed$newcol, split = ":")

    # Create a data frame with each row representing a peak. 
    # "length" determines how many rows per gene depending on how many peaks that gene has
    # "peaks" is the relative location of the peak (relative to what?)
    # Why rename "misc" to "geneWin" here when it's changed back later?
    df.two.peaks = data.frame(chrom = rep(all.two.peaks.bed$chrom, sapply(spikes, length)),
                              start = rep(all.two.peaks.bed$start, sapply(spikes, length)),
                              end = rep(all.two.peaks.bed$end, sapply(spikes, length)),
                              gene = rep(all.two.peaks.bed$gene, sapply(spikes, length)),
                              geneWin = rep(all.two.peaks.bed$misc, sapply(spikes, length)),
                              strand = rep(all.two.peaks.bed$strand, sapply(spikes, length)),
                              peaks = unlist(spikes))
    
    # Change start/end to reflect peak position (strand doesn't matter?)
    df.two.peaks$end <- df.two.peaks$start + as.numeric(as.character(df.two.peaks$peaks))
    df.two.peaks$start <- df.two.peaks$end - 1
    
    # Do up/downstream density and peak height calculations
    upstream.x = upstream.bed(df.two.peaks, denWin)
    up.density = bed6.region.bpQuery.bigWig(bwPlus, bwMinus, upstream.x, op = "avg", abs.value = TRUE, gap.value = 0, bwMap = NULL)
    
    downstream.x = downstream.bed(df.two.peaks, denWin)
    down.density = bed6.region.bpQuery.bigWig(bwPlus, bwMinus, downstream.x, op = "avg", abs.value = TRUE, gap.value = 0, bwMap = NULL)
    
    peak.height = bed6.region.bpQuery.bigWig(bwPlus, bwMinus, df.two.peaks, op = "sum", abs.value = TRUE, gap.value = 0, bwMap = NULL)

    if (densityFilter == TRUE) {
        #next filter based on the density. if down>up, add peak heights and density info for the row.
        df.two.peaks$height[down.density > up.density] = peak.height[down.density > up.density]
        df.two.peaks$up[down.density > up.density] = up.density[down.density > up.density]
        df.two.peaks$down[down.density > up.density] = down.density[down.density > up.density]

        # get all peaks that meet down>up density criteria. this step doesn't keep the "peaks" column
        df.two.peaks = df.two.peaks[!is.na(df.two.peaks$height),c(1,2,3,4,5,6,8,9,10)]
    
        # get predominant TSS by max height (automatically picks first result if multiple max heights, which is why sorting by start is important)
        # df.predominant should be 1 gene = 1 TSS
        temp = df.two.peaks
        df.predominant = do.call(rbind, lapply(split(temp, as.factor(temp$gene)), function(x) {return(x[which.max(x$height),])}))

        # find genes that don't have down>up density TSSs
        two.peaks.default = two.peaks.bed[!(unique(as.character(all.two.peaks.bed$gene)) %in% unique(as.character(df.predominant[,4]))),]
    
        # default to upstream most boundary of exon 1 (based on ref annotation) if the gene is not present
        # the lines won't throw errors even if two.peaks.default has no rows
        two.peaks.default[,2][two.peaks.default$strand == '+'] = two.peaks.default[,2][two.peaks.default$strand == '+'] + tssWin
        two.peaks.default[,3][two.peaks.default$strand == '+'] = two.peaks.default[,2][two.peaks.default$strand == '+'] + 1

        two.peaks.default[,3][two.peaks.default$strand == '-'] = two.peaks.default[,3][two.peaks.default$strand == '-'] - tssWin
        two.peaks.default[,2][two.peaks.default$strand == '-'] = two.peaks.default[,3][two.peaks.default$strand == '-'] - 1
        # check if all columns are present

        two.peaks.default$height <- NA
        two.peaks.default$up <- NA
        two.peaks.default$down <- NA

        colnames(df.predominant) = colnames(two.peaks.default)
    #   print('df.two.peaks')
    #   print(head(df.two.peaks))
        
        df.two.peaks <- rbind(df.predominant, two.peaks.default)
    #   print('two.peaks.default')
    #   print(head(two.peaks.default))
    } else {
        # report height and up/down density for every peak
        df.two.peaks$height = peak.height
        df.two.peaks$up = up.density
        df.two.peaks$down = down.density
        # remove "peaks" column
        df.two.peaks = df.two.peaks[,c(1,2,3,4,5,6,8,9,10)]
        colnames(df.two.peaks) = c("chrom","start","end","gene","misc","strand", "height", "up", "down")
    }

    #rbind the 0, 1, and 2+ dfs
    all.potential.tss <- rbind(no.peaks.bed, one.peak.bed, df.two.peaks)
    return(all.potential.tss)
}

# Other functions to source:
offset.bed <- function(bed) {
    N = dim(bed)[1]
    offby1 = vector(mode="integer", length=N)
    foreach.bed(bed, function(i, chrom, start, end, strand) {
        if (start == end) {
            if (strand == '-') {
                offby1[i] <<- end - 1
            } else {
                offby1[i] <<- start + 1
            } 
        }  
    })
    return(offby1)
}

parse.bed.exon1 <- function(bed.file) {
    bed = read.table(bed.file,
                     col.names=c('chrom', 'start', 'end', 'gene', 'v5', 'strand'))
#identify all strand/gene inconsistencies
    strand.inconsistencies = sapply(strsplit(unique(paste(bed$gene, bed$strand, sep = ':')), split = ":"),
                                    '[', 1)[duplicated(sapply(strsplit(unique(paste(bed$gene, bed$strand, sep = ':')), split = ":"), '[', 1))]
    
#identify all chrY chromosome inconsistencies
    chr.inconsistencies = sapply(strsplit(unique(paste(bed$gene, bed$chr, sep = ':')), split = ":"),
                                 '[', 1)[duplicated(sapply(strsplit(unique(paste(bed$gene, bed$chr, sep = ':')), split = ":"), '[', 1))]
#remove inconsistencies
    bed = bed[!(bed$gene %in%  strand.inconsistencies),]
    bed = bed[(!(bed$gene %in%  chr.inconsistencies) & bed$chr != 'chrY'),]

#identify the upstream-most and downstream-most TSS 
    max.tss.plus = tapply(bed$start, bed$gene, max)[!is.na(tapply(bed$start, bed$gene, max))]
    min.tss.plus = tapply(bed$start, bed$gene, min)[!is.na(tapply(bed$start, bed$gene, min))]

    max.tss.minus = tapply(bed$end, bed$gene, max)[!is.na(tapply(bed$end, bed$gene, max))]
    min.tss.minus = tapply(bed$end, bed$gene, min)[!is.na(tapply(bed$end, bed$gene, min))]

#annotate exon1 range for each gene
    gene = unique(bed$gene)[order(unique(bed$gene))]

#order the strand and chromosome information 
    strand = bed$strand[!duplicated(bed$gene)][order(unique(bed$gene))]
    chrom = bed$chrom[!duplicated(bed$gene)][order(unique(bed$gene))]

#make the dataframe with the correct gene, chrom, and strand info
    uniq.bed = data.frame(chrom, NA, NA, gene, gene, strand, stringsAsFactors = FALSE)
    colnames(uniq.bed) = c('chrom', 'start', 'end', 'gene', 'misc', 'strand')

#annotate the most upstream and most downstream TSS in the second and third columns in a strand-specific manner
    uniq.bed$start[uniq.bed$strand == '+'] <- min.tss.plus[uniq.bed$strand == '+']
    uniq.bed$end[uniq.bed$strand == '+'] <- max.tss.plus[uniq.bed$strand == '+']

    uniq.bed$start[uniq.bed$strand == '-'] <- min.tss.minus[uniq.bed$strand == '-']
    uniq.bed$end[uniq.bed$strand == '-'] <- max.tss.minus[uniq.bed$strand == '-']

    new.start.offset = offset.bed(uniq.bed)

    uniq.bed$end[uniq.bed$strand == '+' & new.start.offset != 0] <- new.start.offset[uniq.bed$strand == '+' & new.start.offset != 0]
    uniq.bed$start[uniq.bed$strand == '-' & new.start.offset != 0] <- new.start.offset[uniq.bed$strand == '-' & new.start.offset != 0]
    return(uniq.bed)
}

top.num <- function(vec.values, top.num.peaks = 20, low.limit.tss.counts = 3) {
# Returns index (position) of top 20 peaks for each bed entry
    subset.len = length(vec.values) - top.num.peaks
    sort.sub = sort(vec.values, partial=subset.len)[subset.len]
    top.20.index = which(vec.values > sort.sub)
    top.20.index = top.20.index[vec.values[top.20.index] >= low.limit.tss.counts]
    return(top.20.index)
}
