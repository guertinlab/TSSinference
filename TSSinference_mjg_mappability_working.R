setwd('~/Desktop/tssInference-master')
library(bigWig)

    
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
    subset.len = length(vec.values) - top.num.peaks
    sort.sub = sort(vec.values, partial=subset.len)[subset.len]
    top.20.index = which(vec.values > sort.sub)
    top.20.index = top.20.index[vec.values[top.20.index] >= low.limit.tss.counts]
    return(top.20.index)
}


#to provide a buffer
#use 5 and 3 prime or up down bed functions.




TSSinference <- function(bed, bw.plus, bw.minus, tssWin = 100, top.num.peaks = 20, low.limit.tss.counts = 3, denWin = 100) {
#buffer    
    bed[,2] = bed[,2] - tssWin
    bed[,3] = bed[,3] + tssWin
#bigwigs
    bwPlus = load.bigWig(bw.plus)
    bwMinus = load.bigWig(bw.minus)

    x.beds = bed6.step.bpQuery.bigWig(bwPlus, bwMinus, bed, step = 1, with.attributes = TRUE)


    v = lapply(x.beds, top.num, top.num.peaks = top.num.peaks, low.limit.tss.counts = low.limit.tss.counts)
    w = as.data.frame(t(sapply(v, "[", i = seq_len(top.num.peaks))))

                                        #now make four data frames: all NA, top.num.peaks - 1 NAs, top.num.peaks -2 NAs, all others
                                        #this is the defaulting to the upstream most TSS instances

                                        #THE DEFAULTING TO UPSTREAM MOST TSS IS NOT CHANGING THE APPROPRIATE END COORDIANTE

    no.peaks = x.beds[is.na(w[,1])]
    no.peaks.bed = bed[is.na(w[,1]),]

    no.peaks.bed[,2][no.peaks.bed$strand == '+'] = no.peaks.bed[,2][no.peaks.bed$strand == '+'] + tssWin
    no.peaks.bed[,3][no.peaks.bed$strand == '-'] = no.peaks.bed[,3][no.peaks.bed$strand == '-'] - tssWin

    no.peaks.bed[,3][no.peaks.bed$strand == '+'] =  no.peaks.bed[,2][no.peaks.bed$strand == '+'] + 1
    no.peaks.bed[,2][no.peaks.bed$strand == '-'] =  no.peaks.bed[,3][no.peaks.bed$strand == '-'] - 1 
    
                                        #done, minus the NA informations
    no.peaks.bed$height <- NA
    no.peaks.bed$up <- NA
    no.peaks.bed$down <- NA
#    print('no.peaks.bed')
#    print(head(no.peaks.bed))
                                        #one peak: if downstream density is > upstream, then use this TSS. otehrwise default to upstream most
    one.peak = w[!is.na(w[,1]) & is.na(w[,2]),]
    one.peak.bed = bed[!is.na(w[,1]) & is.na(w[,2]),]
                                        #dangerous:
    one.peak.bed.transform = one.peak.bed
#this is making a new vector with all the peak coordinates    
    one.peak.bed.transform[,3] = one.peak.bed[,2] + one.peak[,1]
    one.peak.bed.transform[,2] = one.peak.bed.transform[,3] - 1
    
    upstream.x = upstream.bed(one.peak.bed.transform, denWin)

    up.density = bed6.region.bpQuery.bigWig(bwPlus, bwMinus, upstream.x,
                                            op = "avg", abs.value = FALSE, gap.value = 0, bwMap = NULL)
    
    downstream.x = downstream.bed(one.peak.bed.transform, denWin)
    
    down.density = bed6.region.bpQuery.bigWig(bwPlus, bwMinus, downstream.x,
                                              op = "avg", abs.value = FALSE, gap.value = 0, bwMap = NULL)
    
    peak.height= bed6.region.bpQuery.bigWig(bwPlus, bwMinus, one.peak.bed.transform,
                                            op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL)

#use this TRUE FALSE vector and ifelse to default to upstream most or assign a new TSS
#first make them all the upstream most tss, then update it
    
    one.peak.bed[,2][one.peak.bed$strand == '+'] = one.peak.bed[,2][one.peak.bed$strand == '+'] + tssWin
    one.peak.bed[,3][one.peak.bed$strand == '-'] = one.peak.bed[,3][one.peak.bed$strand == '-'] - tssWin
    #one.peak.bed[,2] = one.peak.bed[,2] + tssWin
                                        #one.peak.bed[,3] = one.peak.bed[,3] - tssWin
    
    one.peak.bed[,2][down.density > up.density & one.peak.bed$strand == '+'] <- one.peak.bed.transform[,2][down.density > up.density & one.peak.bed$strand == '+']
    one.peak.bed[,3][down.density > up.density & one.peak.bed$strand == '-'] <- one.peak.bed.transform[,3][down.density > up.density & one.peak.bed$strand == '-']

    one.peak.bed[,3][one.peak.bed$strand == '+'] = one.peak.bed[,2][one.peak.bed$strand == '+'] + 1
    one.peak.bed[,2][one.peak.bed$strand == '-'] = one.peak.bed[,3][one.peak.bed$strand == '-'] - 1
        
    one.peak.bed$height[down.density > up.density] = peak.height[down.density > up.density]
    one.peak.bed$up[down.density > up.density] = up.density[down.density > up.density]
    one.peak.bed$down[down.density > up.density] = down.density[down.density > up.density]
    
#    print('one.peak.bed')
#    print(head(one.peak.bed, 20))
    two.peaks = w[!is.na(w[,3]),]
    two.peaks.bed = bed[!is.na(w[,3]),]
    
    all.two.peaks.bed = cbind(two.peaks.bed, two.peaks)
    
                                        #this is complicated, but quick
    all.two.peaks.bed$newcol <- apply(all.two.peaks.bed[, c(7:(6+top.num.peaks))], 1,
                      function(i){ paste(na.omit(i), collapse = ":") })
    
    spikes <- strsplit(all.two.peaks.bed$newcol, split = ":")
    
    df.two.peaks = data.frame(chrom = rep(all.two.peaks.bed$chrom, sapply(spikes, length)),
                              start = rep(all.two.peaks.bed$start, sapply(spikes, length)),
                              end = rep(all.two.peaks.bed$end, sapply(spikes, length)),
                              gene = rep(all.two.peaks.bed$gene, sapply(spikes, length)),
                              geneWin = rep(all.two.peaks.bed$misc, sapply(spikes, length)),
                              strand = rep(all.two.peaks.bed$strand, sapply(spikes, length)),
                              peaks = unlist(spikes))
    
    df.two.peaks$end <- df.two.peaks$start + as.numeric(as.character(df.two.peaks$peaks))
    df.two.peaks$start <- df.two.peaks$end - 1
    
    upstream.x = upstream.bed(df.two.peaks, denWin)
    
    up.density = bed6.region.bpQuery.bigWig(bwPlus, bwMinus, upstream.x,
                                        op = "avg", abs.value = FALSE, gap.value = 0, bwMap = NULL)
    
    downstream.x = downstream.bed(df.two.peaks, denWin)
    
    down.density = bed6.region.bpQuery.bigWig(bwPlus, bwMinus, downstream.x,
                                          op = "avg", abs.value = FALSE, gap.value = 0, bwMap = NULL)
    
    peak.height= bed6.region.bpQuery.bigWig(bwPlus, bwMinus, df.two.peaks,
                                        op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL)

                                        #next filter based on the density
                                        #find all genes that are present in the original data frame, but not present in this 
    df.two.peaks$height[down.density > up.density] = peak.height[down.density > up.density]
    df.two.peaks$up[down.density > up.density] = up.density[down.density > up.density]
    df.two.peaks$down[down.density > up.density] = down.density[down.density > up.density]
    
                                        #confirm that this works. above deals with care of all none and one instances. two and >2 need to be coded up

    df.two.peaks = df.two.peaks[!is.na(df.two.peaks$height),c(1,2,3,4,5,6,8,9,10)]

    two.peaks.default = two.peaks.bed[!(unique(as.character(all.two.peaks.bed$gene)) %in% unique(as.character(df.two.peaks[,4]))),]
    
                                        #default to upstream most if the gene is not present:
    two.peaks.default[,2][two.peaks.default$strand == '+'] = two.peaks.default[,2][two.peaks.default$strand == '+'] + tssWin
    two.peaks.default[,3][two.peaks.default$strand == '+'] = two.peaks.default[,2][two.peaks.default$strand == '+'] + 1

    two.peaks.default[,3][two.peaks.default$strand == '-'] = two.peaks.default[,3][two.peaks.default$strand == '-'] - tssWin
    two.peaks.default[,2][two.peaks.default$strand == '-'] = two.peaks.default[,3][two.peaks.default$strand == '-'] - 1

    
    two.peaks.default$height <- NA
    two.peaks.default$up <- NA
    two.peaks.default$down <- NA
    
    colnames(df.two.peaks) = colnames(two.peaks.default)
#   print('df.two.peaks')
#    print(head(df.two.peaks))
    df.two.peaks <- rbind(df.two.peaks, two.peaks.default)
#    print('two.peaks.default')
#    print(head(two.peaks.default))

                                        #rbind the 0, 1, and 2+ dfs
#any reason to reorder?
    all.potential.tss <- rbind(no.peaks.bed, one.peak.bed, df.two.peaks)
    return(all.potential.tss)
}

gene.exon1 = parse.bed.exon1('gencode.hg38.firstExon.bed')

x = Sys.time()
potential.tss = TSSinference(gene.exon1, 'HEK293T_dTAG13_PE2_combined_pro_plus.bigWig', 'HEK293T_dTAG13_PE2_combined_pro_minus.bigWig')
y = Sys.time()


x = potential.tss[potential.tss$chrom == 'chr2' & !is.na(potential.tss$height),]

dbscan(x[4000:4008,], eps = 40, minPts=2)



pp = by(x[,c('start', 'height')], as.character(x$gene), dbscan, eps = 40, minPts =2)



dbscan.mods <- function(x) {
    y <- dbscan(x, eps = 40, minPts=2)
    return(y$cluster)
}


pp = by(x[,c('start', 'height')], as.character(x$gene), dbscan.mods)

lp = cbind(x, unlist(pp))



























    
                                        #find clusters of TSS and 
    x.tss = all.potential.tss[!is.na(all.potential.tss$height),]
    
    
    

    
    
    
    
    
    
    
    
    


    






get.peak.from.index <- function(x, new.df) {
    non.na.cols = x[7:length(x[!is.na(x)])]
    print(min(which(is.na(df.two.peaks[,1]))))
    
    new.df[min(which(is.na(df.two.peaks[,1]))),1:6] <- c(x[i,1], as.numeric(as.character(x[i,2])) + as.numeric(as.character(non.na.cols[1])), x[i,3], x[i,4], x[i,5], x[i,6])
    return(new.df)
}

apply(all.two.peaks.bed[1:3,], 1, get.peak.from.index, new.df = df.two.peaks) 

foreach.bed(two.peaks.bed, get.peak.from.index, index.df = two.peaks)

apply(two.peaks, 1, get.peak.from.index) 

    
one.peak.bed.transform[,3] = one.peak.bed[,2] + one.peak[,1]

                                        #the number of rows of a new DF
df <- data.frame(matrix(ncol = 9, nrow = sum(!is.na(two.peaks))))
#need to write a clever function to make a list of data frames or a dataframe output




sum(!is.na(two.peaks))
far.apart = abs(two.peaks[,2] - two.peaks[,1]) >= denWin
                                        #if far apart, calculate density and update either the upstream most TSS or the new tss

#otherwise take the distance between and make a new data frame with the for.each bed function usign the i index to sepcify the upstream and downstream values! look at andres upstream and downstream bed functions



#make a big df with the beginning start and end coordinates and a new column for the index in?
more.peaks = x.beds[!is.na(w[,3])]



get.start <- function(lst) {
    x = attributes(lst)
    return(x$start)
}
get.end <- function(lst) {
    x = attributes(lst)
    return(x$end)
}

get.peaks <- function(in.bed, lst) {
    starts = as.data.frame(sapply(lapply(lst, get.start), "["))
    ends = as.data.frame(sapply(lapply(lst, get.end), "["))
    in.bed$start
}
a = lapply(one.peak, get.start)
b = as.data.frame(sapply(a, "["))
one.peak.bed$peak = b

density.up.down <- function(value.list, bed.input, top.index, INDX, half.window){
    upDen = sum(step.bpQuery.bigWig(bw, chr.val, attr(vector.signal, 'start') + top.index[INDX] - half.window,
                                    attr(vector.signal, 'start') + top.index[INDX] - 1,
                                    step = 1, strand = strand.val, with.attributes = TRUE))/half.window
    downDen = sum(step.bpQuery.bigWig(bw, chr.val, attr(vector.signal, 'start') + top.index[INDX],
                                      attr(vector.signal, 'start') + top.index[INDX] + half.window,
                                      step = 1, strand = strand.val, with.attributes = TRUE))/half.window
}


density.up.down.test <- function(lst) {
    x = attributes(lst)
    starts = as.data.frame(sapply(x, "["))
    
    upDen = sum(step.bpQuery.bigWig(bw, chr.val, attr(vector.signal, 'start') + top.index[INDX] - half.window,
                                    attr(vector.signal, 'start') + top.index[INDX] - 1,
                                    step = 1, strand = strand.val, with.attributes = TRUE))/half.window
    downDen = sum(step.bpQuery.bigWig(bw, chr.val, attr(vector.signal, 'start') + top.index[INDX],
                                      attr(vector.signal, 'start') + top.index[INDX] + half.window,
                                      step = 1, strand = strand.val, with.attributes = TRUE))/half.window
}
                                
a = lapply(one.peak, density.up.down.test)

b = as.data.frame(sapply(a, "["))


                den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, denWin, bw.map = bw.mappability, read.len = read.len)
                if (den[[1]] < den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }

                





#populate each index in the list with the corresponding gene's data frame
list.genes = vector(mode = "list", length = length(x.beds))


as.data.frame(matrix(unlist(listHolder), nrow=length(unlist(listHolder[1])

                                                     

#bed and w have the same number of rows, so indices can apply to one another



#x.beds is like the vec.values
#instead of length(top.20.index) count the NAs.

#need to find clustered peaks and dynamically calcualte densities depedningon how close the peaks are to one another

peaks.bed <- function(bed, bw.plus, bw.minus, top.num.peaks.1 = 20, low.limit.tss.counts.1 = 1, tssWin.1 = 100) {
    bwPlus.1 = load.bigWig(bw.plus)
    bwMinus.1 = load.bigWig(bw.minus)

    list.genes = vector(mode = "list", length = nrow(bed))
                                        #print(list.genes)
#transform bed here with tssWin
#apply bed6.step.bpQuery.bigWig
                                        #apply the sorting, etc to each vector in the list
    lapply(bed6.out, 
    as.data.frame(t(sapply(list.genes, "[", i = seq_len(top.num.peaks.1))))
    foreach.bed(bed, function(i, chrom, start, end, strand, bwPlus = bwPlus.1, bwMinus = bwMinus.1, tssWin = tssWin.1, top.num.peaks = top.num.peaks.1, low.limit.tss.counts = low.limit.tss.counts.1) {
#        substitute:bed6.step.bpQuery.bigWig ?
        if (strand == '+') {
            vec.values = step.bpQuery.bigWig(bwPlus, chrom, start - tssWin,
                                             end + tssWin, step = 1, strand = strand, with.attributes = FALSE)
        } else {
            vec.values = step.bpQuery.bigWig(bwMinus, chrom, start - tssWin,
                                             end + tssWin, step = 1, strand = strand, with.attributes = FALSE)
        }
                                        #print(vec.values)
        subset.len = length(vec.values) - top.num.peaks
        sort.sub = sort(vec.values, partial=subset.len)[subset.len]
        top.20.index = which(vec.values > sort.sub)
        top.20.index = top.20.index[vec.values[top.20.index] >= low.limit.tss.counts]
 #       list.genes[[i]] <<- rep(NA, top.num.peaks)
        list.genes[[i]] <<- top.20.index
        #list.genes[[i]] <<- vec.values
    })
    return(as.data.frame(t(sapply(list.genes, "[", i = seq_len(top.num.peaks.1)))))
}


x = peaks.bed(bed, bw.plus = 'HEK293T_dTAG13_PE2_combined_pro_plus.bigWig',  bw.minus = 'HEK293T_dTAG13_PE2_combined_pro_minus.bigWig')




TSS.infer <- function(bw.plus, bw.minus, bw.map = NULL, read.len = 38, bed.file, tssWin = 100, denWin = 200, top.num.peaks = 20, clustered.peak.distance = 5, low.limit.tss.counts = 3) {
    bed=read.table(bed.file,
                   col.names=c('chrom', 'start', 'end', 'gene', 'v5', 'strand'))
    bwMinus=load.bigWig(bw.minus)
    bwPlus=load.bigWig(bw.plus)
#this is from Luther, I think the biggest thing that it does is address genes that have ambiguous strand annotations:
#need to make it more succinct
    agDf=activeGene(bwPlus,bwMinus,bed,tssWin=tssWin)

#the mappability is not yet tested and may be able to be dropped
    if (!is.null(bw.map)) {
        bw.mappability = load.bwMap(bw.map, read.len = read.len, read.left.edge=TRUE, threshold.fraction = 0)
    } else {
        bw.mappability = NULL
    }
#make a list of dataframes that I will want to combine                                       
    list.genes = vector(mode = "list", length = nrow(bed))
    
    count = 0
    df.fill.out = data.frame(matrix(ncol = 7, nrow = 0))
    denWin = denWin + 1 #this to get a window that is consistent with the user's intention without the need to change density.up.down
    clustered.peak.distance = clustered.peak.distance +1 #this to get a window that is consistent with the user's intention without the need to change density.up.down
    
    for (i in 1:nrow(agDf)) {
        df.fill = data.frame(matrix(ncol = 7, nrow = 0))
        colnames(df.fill) = c('gene', 'chr', 'peak', 'strand', 'up', 'down', 'height') 
        count.1 = 0
 #plus strand
        
        if (agDf$strand[i] == '+') {
            strand.value = '+'
            vec.values = step.bpQuery.bigWig(bwPlus, agDf$chrom[i], agDf$start[i] - tssWin,
                                             agDf$end[i] +tssWin, step = 1, strand = '+', with.attributes = TRUE)
            len.vec.values = length(vec.values)
            subset.len = len.vec.values - top.num.peaks
            sort.sub = sort(vec.values, partial=subset.len)[subset.len]
            top.20.index = which(vec.values > sort.sub)
            top.20.index = top.20.index[vec.values[top.20.index] >= low.limit.tss.counts]
            chr.value = attr(vec.values, 'chrom')
#should we just make a dataframe that is top.num.peaks wide and nrow(bed) long filled with NA? then insert the top.20.index in the respective columns?
#if we do this, then we can operate across the dataframe
#deal with list lengths of 1 differently because there are no neighbors
            if (length(top.20.index) == 1) {
                den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, denWin, bw.map = bw.mappability, read.len = read.len)
                if (den[[1]] < den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
            }#do I just skip the gene if the density downstream is not greater than the upstream?
            
                                        #length of list == 2        
            else if (length(top.20.index) == 2) {
                                        #test case DNAH1, works as intended, but default parameters do not pick up a true TSS            
                if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
                    for (j in 1:2) {
                                        #mappability not yet incorporated
                        den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, denWin, bw.map = bw.mappability, read.len = read.len)
                        if (den[[1]] < den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
                }
                else {           
                    newWin = abs(top.20.index[2] - top.20.index[1])
                    if (newWin > clustered.peak.distance) {
                        for (j in 1:2) {
                        den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, newWin, bw.map = bw.mappability, read.len = read.len)
                            if (den[[1]] < den[[2]]) {
                                count.1 = count.1 + 1
                                df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                            }
                        }
                    } else {
                        for (j in 1:2) {
                            den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance, bw.map = bw.mappability, read.len = read.len)
                            count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
                }
            }
                                        #if vector greater than two the first instance needs special attention because there is only one neighbor            
            else if (length(top.20.index) > 2) {
                if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, denWin, bw.map = bw.mappability, read.len = read.len)
                    if (den[[1]] < den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                    }
                } else {
                    newWin = abs(top.20.index[2] - top.20.index[1])
                    if (newWin > clustered.peak.distance) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, newWin, bw.map = bw.mappability, read.len = read.len)
                        if (den[[1]] < den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                        }
                    }
                else {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, clustered.peak.distance, bw.map = bw.mappability, read.len = read.len)
                    count.1 = count.1 + 1
                    df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
                }
                                        #all internal indicies have two neighbors
                for (j in 2:(length(top.20.index)-1)) {
                    if (abs(top.20.index[j] - top.20.index[j-1])  >= denWin & abs(top.20.index[j] - top.20.index[j+1]) >= denWin) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, denWin, bw.map = bw.mappability, read.len = read.len)
                        if (den[[1]] < den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
                    else {
                        newWin = min(c(abs(top.20.index[j-1] - top.20.index[j]), abs(top.20.index[j+1] - top.20.index[j])))
                        if (newWin > clustered.peak.distance) {
                            den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, newWin, bw.map = bw.mappability, read.len = read.len)
                            if (den[[1]] < den[[2]]) {
                                count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                            }
                        } else {
                            count.1 = count.1 + 1
                            den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance, bw.map = bw.mappability, read.len = read.len)
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
            }
                                        #THE LAST ENTRY IN THE VECTOR ONLY HAS ONE NEIGHBOR, SO TREAT THIS SPECIAL CASE
                if (abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)]) >= denWin) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), denWin, bw.map = bw.mappability, read.len = read.len)
                    if (den[[1]] < den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                    }
                }
            else {
                newWin = abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)])
                if (newWin > clustered.peak.distance) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), newWin, bw.map = bw.mappability, read.len = read.len)
                    if (den[[1]] < den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values)
                    }
                } else {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), clustered.peak.distance, bw.map = bw.mappability, read.len = read.len)
                    count.1 = count.1 + 1
                    df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values)
                }
            }     
            }
        }
        else if (agDf$strand[i] == '-') {
            strand.value = '-'
            vec.values = step.bpQuery.bigWig(bwMinus, agDf$chrom[i], agDf$start[i] - tssWin,
                                             agDf$end[i] +tssWin, step = 1, strand = '-', with.attributes = TRUE)
            len.vec.values = length(vec.values)
            subset.len = len.vec.values - top.num.peaks
            sort.sub = sort(vec.values, partial=subset.len)[subset.len]
            top.20.index = which(vec.values > sort.sub)
            top.20.index = top.20.index[vec.values[top.20.index] >= low.limit.tss.counts]
            chr.value = attr(vec.values, 'chrom')
            
            if (length(top.20.index) == 1) {
                den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, denWin, bw.map = bw.mappability, read.len = read.len)
                                        #simply swap the up and down logical operator (less than to greater than)
                if (den[[1]] > den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)  
                }#what do we do if the density is not greater downstream.
            }
            else if (length(top.20.index) == 2) {
                if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
                    for (j in 1:2) {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, denWin, bw.map = bw.mappability, read.len = read.len)
                        if (den[[1]] > den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
                }
                else {
                newWin = abs(top.20.index[2] - top.20.index[1])
                if (newWin > clustered.peak.distance) {
                    for (j in 1:2) {
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, newWin, bw.map = bw.mappability, read.len = read.len)
                        if (den[[1]] > den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
                } else {
                    for (j in 1:2) {
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance, bw.map = bw.mappability, read.len = read.len)
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                    }
                }
                }
            }
        else if (length(top.20.index) > 2) {
            if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
                den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, denWin, bw.map = bw.mappability, read.len = read.len)
                if (den[[1]] > den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
            } else {
                newWin = abs(top.20.index[2] - top.20.index[1])
                if (newWin > clustered.peak.distance) {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, newWin, bw.map = bw.mappability, read.len = read.len)
                    if (den[[1]] > den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                    }
                } else {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, clustered.peak.distance, bw.map = bw.mappability, read.len = read.len)
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
            }
            for (j in 2:(length(top.20.index)-1)) {
                if (abs(top.20.index[j] - top.20.index[j-1])  >= denWin & abs(top.20.index[j] - top.20.index[j+1]) >= denWin) {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, denWin, bw.map = bw.mappability, read.len = read.len)
                    if (den[[1]] > den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                    }
                }
                else {
                    newWin = min(c(abs(top.20.index[j-1] - top.20.index[j]), abs(top.20.index[j+1] - top.20.index[j])))
                    if (newWin > clustered.peak.distance) {
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, newWin, bw.map = bw.mappability, read.len = read.len)
                        if (den[[1]] > den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    } else {
                        count.1 = count.1 + 1
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance, bw.map = bw.mappability, read.len = read.len)
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values) 
                    }
                }
            }
                                        #THE LAST ENTRY IN THE VECTOR ONLY HAS ONE NEIGHBOR, SO TREAT THIS SPECIAL CASE
            if (abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)]) >= denWin) {
                den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), denWin, bw.map = bw.mappability, read.len = read.len)
                if (den[[1]] > den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                }
            }
            else {
                newWin = abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)])
                if (newWin > clustered.peak.distance) {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), newWin, bw.map = bw.mappability, read.len = read.len)
                    if (den[[1]] > den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                    }
                } else {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), clustered.peak.distance, bw.map = bw.mappability, read.len = read.len)
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                }
            }     
        }
        }
                                        #this deals with issues where there is no entry for a gene, defaults to upstream most TSS    
        if (nrow(df.fill) == 0 & strand.value == '+') {
            df.fill[1, 3] = attr(vec.values, 'start') + tssWin
            df.fill[,1] = agDf[i, 1]
            df.fill[,2] = attr(vec.values, 'chrom')
            df.fill[,4] = agDf[i, 5]  
            df.fill[1, 5] = NA
            df.fill[1, 6] = NA
            df.fill[1, 7] = NA
        }
        else if (nrow(df.fill) == 0 & strand.value == '-') {                               
            df.fill[1, 3] = attr(vec.values, 'end') - tssWin
            df.fill[,1] = agDf[i, 1]
            df.fill[,2] = attr(vec.values, 'chrom')
            df.fill[,4] = agDf[i, 5]  
            df.fill[1, 5] = NA
            df.fill[1, 6] = NA
            df.fill[1, 7] = NA
        }
        df.fill.out = rbind(df.fill.out, df.fill)
    }
    if (!is.null(bw.mappability)) {
        df.fill.out$map = ifelse((df.fill.out$up == 0 & df.fill.out$down == 0 & !is.na(df.fill.out$up)), "not mappable", "mappable")
    }
                                        #filter df.fill.out to provide a single predominant TSS per gene
    temp = df.fill.out
    temp$height[is.na(temp$height)] <- -2
    df.predominant = do.call(rbind, lapply(split(temp, as.factor(temp$gene)), function(x) {return(x[which.max(x$height),])}))
    df.predominant$height[df.predominant$height == -2] <- NA
    return(list(df.fill.out, df.predominant))
}

                                        #I had labeled the input files incorrectly, plus is minus and minus is plus
bw.mappability.1 = load.bwMap('hg38_38mers.bigWig', read.len = 38, read.left.edge=TRUE, threshold.fraction = 0)

step.bpQuery.bwMap(bw.mappability.1, 'chr12', 103932303, 103932403,
                   step = 1, op = "thresh", strand = '+', with.attributes = TRUE)



tss.test = TSS.infer(bw.plus = 'H9_minus_PE2_merged.bigWig', bw.minus = 'H9_plus_PE2_merged.bigWig', bw.map = 'hg38_38mers.bigWig', bed.file = 'gencode.hg38.firstExon.bed')

tss.test.no.map = TSS.infer(bw.plus = 'HEK293T_dTAG13_PE2_combined_pro_plus.bigWig', bw.minus = 'HEK293T_dTAG13_PE2_combined_pro_minus.bigWig', bed.file = 'gencode.hg38.firstExon.bed')

save(tss.test, tss.test.no.map, file = 'tssinference.Rdata')












activeGene <- function(bwPlus, bwMinus, bed, tssWin){
  ##--Arg--##
  #bed: bed as a data frame of 1st exons in a gene
  #bwPlus: bigWig file for the plus strand
  #bwMinus: bigWig file for the minus strand
  #tssWin: integer number of bp to look at beforegene start

  #cat('Debugging',sep='')
  # list unique genes in bed file
  gene=unique(bed$gene)
  # initiate empty vectors
  chrom=as.character(c())
  start=as.integer(c())# start position of gene
  end=as.integer(c())# end position of gene
  strand=as.character(c())# strand to use for preak calculation

  plusCnt=as.integer(c())# number of exons on plus strand in gene
  minusCnt=as.integer(c())#  number of exons on minus strand in gene
    # calculated using bigWig functions
  active=as.logical(c())# logical whether gene is active
  region=as.integer(c())# integer of pileups in region

  warning=as.integer(c())# vector of integers describing warning
  warningFlag=as.logical(c())# logical flag for sorting

  # loop through unique genes
  for(i in gene){
    #reset gene warning vector
    gWarning=as.integer(c())

    # create bed of gene isoforms
    gBed=bed[bed$gene==i,]

    ################
    # determine chrom
    gChrom=as.character(unique(gBed$chrom))
    if(length(gChrom)>1){
      chrom=append(chrom,'chrX')# append to chrom vector
      gWarning=append(gWarning,1)
    }
    else if(length(gChrom)==1){
      chrom=append(chrom,as.character(gChrom))
    }

    ################
    # determine strandness
    pC=dim(gBed[gBed$gene==i&gBed$strand=='+',])[1]
    mC=dim(gBed[gBed$gene==i&gBed$strand=='-',])[1]
    if(pC>mC){
      gStrand='+'
    }
    else if(mC>pC){
      gStrand='-'
    }
    #-----------------------#
    else if(mC==pC){
      gbedPlus=gBed[gBed$strand=='+']
      #num of bp
      gbedPlusbp=max(gBedPlus$start)-min(gBedPlus$start)
      gbedMinus=gBed[gBed$strand=='-']
      gbedMinusbp=max(gBedMinus$start)-min(gBedMinus$start)
      if(gbedPlusbp>gbedMinusbp){
        gStrand='+'
      }
      else if(gbedPlusbp<gbedMinusbp){
        gStrand='-'
      }
      else if(gbedPlusbp==gbedMinusbp){
        gStrand=NA
      }
    }# need to test
    strand=append(strand,gStrand)
    #-----------------------#

    #add to vectors
    plusCnt=append(plusCnt,pC)
    minusCnt=append(minusCnt,mC)
    #if(gStrand==NA){gWarnings=append(gWarnings,2}
    ################
    # determine start, end
    gBedStarts=c()
    for(row in 1:nrow(gBed)){
      if(as.character(gBed[row,'strand'])=='+'){
        gBedStarts=append(gBedStarts,as.integer(gBed[row,'start']))
      }
      else if(as.character(gBed[row,'strand'])=='-'){
        gBedStarts=append(gBedStarts,as.integer(gBed[row,'end']))
      }
    }
    # append gene start and end
    gStart=min(gBedStarts)
    start=append(start,min(gBedStarts))
    gEnd=max(gBedStarts)
    end=append(end,max(gBedStarts))

    if(length(gChrom)>1){gChrom=gChrom[1]}
    #use start, end to determine active
    if(gStart!=gEnd){
      if(gStrand=='+'){
        gRegion=region.bpQuery.bigWig(bwPlus,start=gStart-tssWin,
                                  end=gEnd+tssWin,
                                  chrom=gChrom)
                                }
      if(gStrand=='-'){
        gRegion=region.bpQuery.bigWig(bwMinus,start=gStart-tssWin,
                                    end=gEnd+tssWin,
                                    chrom=gChrom)
                                }}
    else if(gStart==gEnd){gRegion=0}
    if(gRegion>0){
        active=append(active,TRUE)
        region=append(region, gRegion)
    }
    else if(gRegion==0){
        active=append(active,FALSE)
        region=append(region, gRegion)
    }

  }
  dfReturn=data.frame(gene,chrom,start,end,strand,
    active,plusCnt, minusCnt, region,stringsAsFactors = FALSE)
  #add evaluate column
  dfReturn['evaluate']=TRUE
  #add geneWin length of bp in gene window
  dfReturn['geneWin']=dfReturn$end-dfReturn$start

  #remove unwanted genes from being evaluated
  dfReturn$evaluate[dfReturn$active==FALSE]=FALSE
  dfReturn$evaluate[(dfReturn$strand!='+'&dfReturn$strand!='-')]=FALSE
  dfReturn$evaluate[dfReturn$geneWin>2500000]=FALSE

  #return dataframe
  return(dfReturn)
}



#rewrite activeGene       


#since start sites are empirically clustered, then it would make sense to favor the highest peak in a cluster of peaks.




bed.test=read.table('gencode.hg38.firstExon.bed',
                   col.names=c('chrom', 'start', 'end', 'gene', 'v5', 'strand'))
bwPlus=load.bigWig('HEK293T_dTAG13_PE2_combined_pro_plus.bigWig')
bwMinus= load.bigWig('HEK293T_dTAG13_PE2_combined_pro_minus.bigWig')

agDf=activeGene(bwPlus ,bwMinus, bed.test, tssWin=tssWin)

                                        #I think this does activeGene much quicker

bed.file = 'gencode.hg38.firstExon.bed'

#quickly adjusts the start and end columns in a strand-aware manner for genes that have a single exon
#uses foreach.bed from bigWig library
