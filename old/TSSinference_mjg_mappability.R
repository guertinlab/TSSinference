setwd('~/Desktop/tssInference-master')
#activeGene reproduced below
#source('https://raw.githubusercontent.com/lutus/tssInference/master/activeGene.R')
library(bigWig)

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

density.up.down <- function(bw, chr.val, vector.signal, strand.val, top.index, INDX, half.window, bw.map = NULL){
    if (is.null(bw.map)) {
        upDen = sum(step.bpQuery.bigWig(bw, chr.val, attr(vector.signal, 'start') + top.index[INDX] - half.window,
                                    attr(vector.signal, 'start') + top.index[INDX] - 1,
                                    step = 1, strand = '+', with.attributes = TRUE))/half.window
        downDen = sum(step.bpQuery.bigWig(bw, chr.val, attr(vector.signal, 'start') + top.index[INDX],
                                      attr(vector.signal, 'start') + top.index[INDX] + half.window,
                                      step = 1, strand = strand.val, with.attributes = TRUE))/half.window

    } else {
#the *3 is a buffer, if three times the window does not result in enough mappable bases, then consider defaulting to upstream most TSS
    test.map.up = step.bpQuery.bigWig(bw.map, chr.val, attr(vector.signal, 'start') + top.index[INDX] - half.window*3,
                                      attr(vector.signal, 'start') + top.index[INDX] -1,
                                      step = 1, with.attributes = TRUE)
    test.map.down = step.bpQuery.bigWig(bw.map, chr.val, attr(vector.signal, 'start') + top.index[INDX],
                                        attr(vector.signal, 'start') + top.index[INDX] + half.window*3 -1,
                                        step = 1, with.attributes = TRUE)
    index.up = match(half.window, cumsum(abs(test.map.up -1)))
    index.down = match(half.window, cumsum(abs(test.map.down -1)))
        if (!is.na(index.up) & !is.na(index.down)) {
            upDen = sum(step.bpQuery.bigWig(bw, chr.val, attr(vector.signal, 'start') + top.index[INDX] - index.up,
                                            attr(vector.signal, 'start') + top.index[INDX] - 1,
                                            step = 1, strand = '+', with.attributes = TRUE))/(half.window-1)
            downDen = sum(step.bpQuery.bigWig(bw, chr.val, attr(vector.signal, 'start') + top.index[INDX],
                                              attr(vector.signal, 'start') + top.index[INDX] + index.down -1,
                                              step = 1, strand = strand.val, with.attributes = TRUE))/(half.window-1)
        } else {
            upDen = 0
            downDen = 0
#            print(chr.val)
#            print(attr(vector.signal, 'start') + top.index[INDX])
        }
    }
    return(list(upDen, downDen))
}

add.to.fill <- function(orig.df, df.filler, dense.u.d, counter, line.num, top.index, IDX, vector.signal) {
    df.filler[counter, 3] = attr(vector.signal, 'start') + top.index[IDX]
    df.filler[,1] = orig.df[line.num, 1]
    df.filler[,2] = attr(vector.signal, 'chrom')
    df.filler[,4] = orig.df[line.num, 5]  
    df.filler[counter, 5] = dense.u.d[[1]]
    df.filler[counter, 6] = dense.u.d[[2]]
    df.filler[counter, 7] = vector.signal[top.index[IDX]]
    return(df.filler)
}

add.to.fill.minus <- function(orig.df, df.filler, dense.u.d, counter, line.num, top.index, IDX, vector.signal) {
    df.filler[counter, 3] = attr(vector.signal, 'end') - (length(vector.signal) - top.index[IDX])
    df.filler[,1] = orig.df[line.num, 1]
    df.filler[,2] = attr(vector.signal, 'chrom')
    df.filler[,4] = orig.df[line.num, 5]  
    df.filler[counter, 5] = dense.u.d[[1]]
    df.filler[counter, 6] = dense.u.d[[2]]
    df.filler[counter, 7] = vector.signal[top.index[IDX]]
    return(df.filler)
}



TSS.infer <- function(bw.plus, bw.minus, bw.map = NULL, bed.file, tssWin = 100, denWin = 200, top.num.peaks = 20, clustered.peak.distance = 5, low.limit.tss.counts = 3) {
#bw.plus = 'H9_plus_PE2_merged.bigWig', bw.minus = 'H9_minus_PE2_merged.bigWig', bw.map = 'hg38_38mers.bigWig', bed.file = 'gencode.hg38.firstExon.bed'
    bed=read.table(bed.file,
                   col.names=c('chrom', 'start', 'end', 'gene', 'v5', 'strand'))
    bwMinus=load.bigWig(bw.minus)
    bwPlus=load.bigWig(bw.plus)

    agDf=activeGene(bwPlus,bwMinus,bed,tssWin=tssWin)

    if (!is.null(bw.map)) {
        bw.mappability = load.bigWig(bw.map)
    } else {
        bw.mappability = NULL
    }
    
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
                                        #deal with list lengths of 1 differently because there are no neighbors
            if (length(top.20.index) == 1) {
                den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, denWin, bw.map = bw.mappability)
                if (den[[1]] < den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
            }
                                        #length of list == 2        
            else if (length(top.20.index) == 2) {
                                        #test case DNAH1, works as intended, but default parameters do not pick up a true TSS            
                if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
                    for (j in 1:2) {
                                        #mappability not yet incorporated
                        den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, denWin, bw.map = bw.mappability)
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
                        den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, newWin, bw.map = bw.mappability)
                            if (den[[1]] < den[[2]]) {
                                count.1 = count.1 + 1
                                df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                            }
                        }
                    } else {
                        for (j in 1:2) {
                            den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance, bw.map = bw.mappability)
                            count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
                }
            }
                                        #if vector greater than two the first instance needs special attention because there is only one neighbor            
            else if (length(top.20.index) > 2) {
                if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, denWin, bw.map = bw.mappability)
                    if (den[[1]] < den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                    }
                } else {
                    newWin = abs(top.20.index[2] - top.20.index[1])
                    if (newWin > clustered.peak.distance) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, newWin, bw.map = bw.mappability)
                        if (den[[1]] < den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                        }
                    }
                else {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, clustered.peak.distance, bw.map = bw.mappability)
                    count.1 = count.1 + 1
                    df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
                }
                                        #all internal indicies have two neighbors
                for (j in 2:(length(top.20.index)-1)) {
                    if (abs(top.20.index[j] - top.20.index[j-1])  >= denWin & abs(top.20.index[j] - top.20.index[j+1]) >= denWin) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, denWin, bw.map = bw.mappability)
                        if (den[[1]] < den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
                    else {
                        newWin = min(c(abs(top.20.index[j-1] - top.20.index[j]), abs(top.20.index[j+1] - top.20.index[j])))
                        if (newWin > clustered.peak.distance) {
                            den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, newWin, bw.map = bw.mappability)
                            if (den[[1]] < den[[2]]) {
                                count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                            }
                        } else {
                            count.1 = count.1 + 1
                            den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance, bw.map = bw.mappability)
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
            }
                                        #THE LAST ENTRY IN THE VECTOR ONLY HAS ONE NEIGHBOR, SO TREAT THIS SPECIAL CASE
                if (abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)]) >= denWin) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), denWin, bw.map = bw.mappability)
                    if (den[[1]] < den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                    }
                }
            else {
                newWin = abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)])
                if (newWin > clustered.peak.distance) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), newWin, bw.map = bw.mappability)
                    if (den[[1]] < den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values)
                    }
                } else {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), clustered.peak.distance, bw.map = bw.mappability)
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
                den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, denWin, bw.map = bw.mappability)
                                        #simply swap the up and down logical operator (less than to greater than)
                if (den[[1]] > den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)  
                }
            }
            else if (length(top.20.index) == 2) {
                if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
                    for (j in 1:2) {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, denWin, bw.map = bw.mappability)
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
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, newWin, bw.map = bw.mappability)
                        if (den[[1]] > den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
                } else {
                    for (j in 1:2) {
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance, bw.map = bw.mappability)
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                    }
                }
                }
            }
        else if (length(top.20.index) > 2) {
            if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
                den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, denWin, bw.map = bw.mappability)
                if (den[[1]] > den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
            } else {
                newWin = abs(top.20.index[2] - top.20.index[1])
                if (newWin > clustered.peak.distance) {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, newWin, bw.map = bw.mappability)
                    if (den[[1]] > den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                    }
                } else {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, clustered.peak.distance, bw.map = bw.mappability)
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
            }
            for (j in 2:(length(top.20.index)-1)) {
                if (abs(top.20.index[j] - top.20.index[j-1])  >= denWin & abs(top.20.index[j] - top.20.index[j+1]) >= denWin) {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, denWin, bw.map = bw.mappability)
                    if (den[[1]] > den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                    }
                }
                else {
                    newWin = min(c(abs(top.20.index[j-1] - top.20.index[j]), abs(top.20.index[j+1] - top.20.index[j])))
                    if (newWin > clustered.peak.distance) {
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, newWin, bw.map = bw.mappability)
                        if (den[[1]] > den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    } else {
                        count.1 = count.1 + 1
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance, bw.map = bw.mappability)
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values) 
                    }
                }
            }
                                        #THE LAST ENTRY IN THE VECTOR ONLY HAS ONE NEIGHBOR, SO TREAT THIS SPECIAL CASE
            if (abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)]) >= denWin) {
                den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), denWin, bw.map = bw.mappability)
                if (den[[1]] > den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                }
            }
            else {
                newWin = abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)])
                if (newWin > clustered.peak.distance) {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), newWin, bw.map = bw.mappability)
                    if (den[[1]] > den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                    }
                } else {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), clustered.peak.distance, bw.map = bw.mappability)
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
tss.test = TSS.infer(bw.plus = 'H9_minus_PE2_merged.bigWig', bw.minus = 'H9_plus_PE2_merged.bigWig', bw.map = 'hg38_38mers.bigWig', bed.file = 'gencode.hg38.firstExon.bed')

tss.test.no.map = TSS.infer(bw.plus = 'H9_plus_PE2_merged.bigWig', bw.minus = 'H9_minus_PE2_merged.bigWig', bed.file = 'gencode.hg38.firstExon.bed')

save(tss.test, tss.test.no.map, file = 'tssinference.Rdata')









