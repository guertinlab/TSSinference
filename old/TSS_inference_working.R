setwd('~/Desktop/tssInference-master')
source('https://raw.githubusercontent.com/lutus/tssInference/master/activeGene.R')

library(bigWig)
loc='gencode.hg38.firstExon.bed.gz'
bed=read.table(gzfile(loc),
             col.names=c('chrom', 'start', 'end', 'gene', 'v5', 'strand'))

bwMinusLoc='H9_plus_PE2_merged.bigWig'
bwPlusLoc='H9_minus_PE2_merged.bigWig'
bwMinus=load.bigWig(bwMinusLoc)
bwPlus=load.bigWig(bwPlusLoc)

agDf=activeGene(bwPlus,bwMinus,bed,tssWin=100)

#This is named tssWin, but it needs to be changed to denWin


#agDf.stored = agDf
                                        #agDf = agDf[agDf$gene == 'APPL2',]
                                        #agDf = agDf[1:500,]

                                        #for each gene I would like the following for the top 20 peaks
                                        #gene chr peak strand up.density down.density height
                                      #then we can weight height and density if we carry these values into a DF

density.up.down <- function(bw, chr.val, vector.signal, strand.val, top.index, INDX, half.window ){
    upDen = sum(step.bpQuery.bigWig(bw, chr.val, attr(vector.signal, 'start') + top.index[INDX] - half.window,
                                    attr(vector.signal, 'start') + top.index[INDX] - 1,
                                    step = 1, strand = '+', with.attributes = TRUE))/(half.window-1)

    downDen = sum(step.bpQuery.bigWig(bw, chr.val, attr(vector.signal, 'start') + top.index[INDX],
                                      attr(vector.signal, 'start') + top.index[INDX] + half.window -1,
                                      step = 1, strand = strand.val, with.attributes = TRUE))/(half.window-1)

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
                                 
t0=Sys.time()


#funtion arguments
tssWin = 100
denWin = 200 #not inclusive in the calculation, so this woudl be 199
top.num.peaks = 20
clustered.peak.distance = 5
low.limit.tss.counts = 3


#would be where the function starts
count = 0
df.fill.out = data.frame(matrix(ncol = 7, nrow = 0))
denWin = denWin + 1 #this to get a window that is consistent with the user's intention without the need to change density.up.down
clustered.peak.distance = clustered.peak.distance +1 #this to get a window that is consistent with the user's intention without hte need to change density.up.down

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
#test case is ALPL. although this is not the correct TSS with default parameters, the statement is working as intended
        if (length(top.20.index) == 1) {
            den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, denWin)
            if (den[[1]] < den[[2]]) {
                count.1 = count.1 + 1
                df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
            }
        }
#something needs to be doen if the densisities are not < and greater than  
#another example like above woule be is the max 
#length of list == 2        
        else if (length(top.20.index) == 2) {
#test case DNAH1, works as intended, but default parameters do not pick up a true TSS            
            if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
                for (j in 1:2) {
#mappability not yet incorporated
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, denWin)
                    if (den[[1]] < den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                    }
                }
            }
            else {
#test case was SPRY1. This is a nice example where that annotated TSS is picked up
#There is an alternative TSS that is not picked up
#the density calculation was including the downstream peak, so I modififed the function.                
                newWin = abs(top.20.index[2] - top.20.index[1])
                if (newWin > clustered.peak.distance) {
                    for (j in 1:2) {
#mappability not yet incorporated
                        den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, newWin)
                        if (den[[1]] < den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
                } else {
#HERE IS A BIG ISSUE. IF THERE ARE SEVERAL PEAKS NEARBY ONE ANOTHER, THEN DENSITY UP and DOWN become MEANINGLESS.
#THE RULES ARE MORE RELEVANT HERE, BUT ALSO NEED TO BE APPLIED TO INDEX VECTORS OF length == 2
#FOR NOW I WANT TO WRITE SOMETHING THAT JUST LEAVES all these CLUSTERED CLOSE PEAKS IN THE DATA FRAME REGARDLESS OF UP AND DOWN DENSITIES
#MAYBE A USER_IDENTIFIED VALUE FOR THIS (clustered.peak.distance). IF THE PEAKS ARE CLOSER THAN 5 BASES from EACH Other, then WE CAN keep THEM all
#AND CALCULATE density in the user-defined clustered.peak.distance
#THIS SPECIAL CASE NEEDS TO BE APPLIED TO OTHER SMALL DISTANCES
#TMEM151A is the test case. this find a few new TSS within short distances from one another. 
                    for (j in 1:2) {
#mappability not yet incorporated
                        den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance)
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                    }
                }
            }
        }
#if vector greater than two the first instance needs special attention because there is only one neighbor            
        else if (length(top.20.index) > 2) {
#FHAD1 is a test case and a novel confident TSS
            if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
#mappability not yet incorporated
                den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, denWin)
                if (den[[1]] < den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
            } else {
                newWin = abs(top.20.index[2] - top.20.index[1])
                if (newWin > clustered.peak.distance) {
# test case  KLHL17
#mappability not yet incorporated
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, newWin)
                    if (den[[1]] < den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                    }
                }
                else {
#PLEKHN1 is a test case and identifies a new confident TSS                    
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, 1, clustered.peak.distance)
                    count.1 = count.1 + 1
                    df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
            }
#all internal indicies have two neighbors
            for (j in 2:(length(top.20.index)-1)) {
                if (abs(top.20.index[j] - top.20.index[j-1])  >= denWin & abs(top.20.index[j] - top.20.index[j+1]) >= denWin) {
#mappability not yet incorporated
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, denWin)
                    if (den[[1]] < den[[2]]) {
#ESPN is a test case
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                    }
                }
                else {
                    newWin = min(c(abs(top.20.index[j-1] - top.20.index[j]), abs(top.20.index[j+1] - top.20.index[j])))
                    if (newWin > clustered.peak.distance) {
                        den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, newWin)
                        if (den[[1]] < den[[2]]) {
  #B3GALT6
                            count.1 = count.1 + 1
                            df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    } else {
#KLHL17 example again
                        count.1 = count.1 + 1
                        den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance)
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                    }
                }
            }
#THE LAST ENTRY IN THE VECTOR ONLY HAS ONE NEIGHBOR, SO TREAT THIS SPECIAL CASE
            if (abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)]) >= denWin) {
#mappability not yet incorporated
                den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), denWin)
                if (den[[1]] < den[[2]]) {
#PRKCZ                     
                    count.1 = count.1 + 1
                    df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                }
            }
            else {
                newWin = abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)])
                if (newWin > clustered.peak.distance) {
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), newWin)
                    if (den[[1]] < den[[2]]) {
#PHF13                        
                        count.1 = count.1 + 1
                        df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values)
                    }
                } else {
#ISG15
                    den = density.up.down(bwPlus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), clustered.peak.distance)
                    count.1 = count.1 + 1
                    df.fill = add.to.fill(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values)
                }
            }     
        }
    }
#STOPPED HERE ON MAY 9, 2020    
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
            den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, denWin)
            #simply swap the up and down logical operator (less than to greater than)
            if (den[[1]] > den[[2]]) {
                count.1 = count.1 + 1
                df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)  
            }
        }
#SOMETHING TO DO IF THE DENSITIES ARE NOT GREATER LESS...        
        else if (length(top.20.index) == 2) {
            if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
                for (j in 1:2) {
                                        #mappability not yet incorporated
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, denWin)
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
#mappability not yet incorporated
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, newWin)
                        if (den[[1]] > den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    }
                } else {
                    for (j in 1:2) {
#mappability not yet incorporated
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance)
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                    }
                }
            }
        }
        else if (length(top.20.index) > 2) {
            if (abs(top.20.index[2] - top.20.index[1]) >= denWin) {
#mappability not yet incorporated
                den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, denWin)
                if (den[[1]] > den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
            } else {
                newWin = abs(top.20.index[2] - top.20.index[1])
                if (newWin > clustered.peak.distance) {
#mappability not yet incorporated
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, newWin)
                    if (den[[1]] > den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                    }
                } else {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, 1, clustered.peak.distance)
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, 1, vec.values)
                }
            }
#all internal indicies have two neighbors
            for (j in 2:(length(top.20.index)-1)) {
                if (abs(top.20.index[j] - top.20.index[j-1])  >= denWin & abs(top.20.index[j] - top.20.index[j+1]) >= denWin) {
#mappability not yet incorporated
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, denWin)
                    if (den[[1]] > den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                    }
                }
                else {
                    newWin = min(c(abs(top.20.index[j-1] - top.20.index[j]), abs(top.20.index[j+1] - top.20.index[j])))
                    if (newWin > clustered.peak.distance) {
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, newWin)
                        if (den[[1]] > den[[2]]) {
                            count.1 = count.1 + 1
                            df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values)
                        }
                    } else {
                        count.1 = count.1 + 1
                        den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, j, clustered.peak.distance)
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, j, vec.values) 
                    }
                }
            }
#THE LAST ENTRY IN THE VECTOR ONLY HAS ONE NEIGHBOR, SO TREAT THIS SPECIAL CASE
            if (abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)]) >= denWin) {
#mappability not yet incorporated
                den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), denWin)
                if (den[[1]] > den[[2]]) {
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                }
            }
            else {
                newWin = abs(top.20.index[length(top.20.index)-1] - top.20.index[length(top.20.index)])
                if (newWin > clustered.peak.distance) {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), newWin)
                    if (den[[1]] > den[[2]]) {
                        count.1 = count.1 + 1
                        df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                    }
                } else {
                    den = density.up.down(bwMinus, chr.value, vec.values, strand.value, top.20.index, length(top.20.index), clustered.peak.distance)
                    count.1 = count.1 + 1
                    df.fill = add.to.fill.minus(agDf, df.fill, den, count.1, i, top.20.index, length(top.20.index), vec.values) 
                }
            }     
        }
    }
#this deals with issues where there is no entry for a gene, defaults to upstream most TSS    
    if (nrow(df.fill) == 0 & strand.value == '+') {
            #print(attr(vec.values, 'start'))                           #if these
            df.fill[1, 3] = attr(vec.values, 'start') + tssWin
            df.fill[,1] = agDf[i, 1]
            df.fill[,2] = attr(vec.values, 'chrom')
            df.fill[,4] = agDf[i, 5]  
            df.fill[1, 5] = NA
            df.fill[1, 6] = NA
            df.fill[1, 7] = NA
    }
    else if (nrow(df.fill) == 0 & strand.value == '-') {                                    #if these
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

                                     


t1=Sys.time()
dt=difftime(t1,t0,units='mins')
dt

df.fill.out[1:100,]

#PRDM16
#APPL2
#stop here



bw.map = load.bigWig('hg38_38mers.bigWig')

test.map.up = step.bpQuery.bigWig(bw.map, chrm, peakBp - littleBuffer, peakBp-1, step = 1, with.attributes = TRUE)
test.map.down = step.bpQuery.bigWig(bw.map, chrm, peakBp + 1, peakBp + littleBuffer, step = 1, with.attributes = TRUE)

index.up = match(denWindow, cumsum(abs(test.map.up -1)))
index.down = match(denWindow, cumsum(abs(test.map.down -1)))

upDen=region.bpQuery.bigWig(bwPlus,chrm, peakBp-index.up, peakBp-1, op='sum')/denWindow

dwnDen=region.bpQuery.bigWig(bwPlus,chrm,peakBp+1,peakBp+index.down, op='sum')/denWindow











