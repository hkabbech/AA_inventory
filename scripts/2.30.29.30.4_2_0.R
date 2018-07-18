setwd(paste(getwd(),'/superfamilies/2.30.29.30.4_2_0/tables_old/',sep=''))
NB1 = 'loop_N_B1_aa.csv'
B1 = 'fold_B1_aa.csv'
B1B2 = 'loop_B2_B3_aa.csv'
B2 = 'fold_B2_aa.csv'
B2B3 = 'loop_B2_B3_aa.csv'
B3 = 'fold_B3_aa.csv'
B3B4 = 'loop_B3_B4_aa.csv'
B4 = 'fold_B4_aa.csv'
B4B5 = 'loop_B4_B5_aa.csv'
B5 = 'fold_B5_aa.csv'
B5B6 = 'loop_B5_B6_aa.csv'
B6 = 'fold_B6_aa.csv'
B6B7 = 'loop_B6_B7_aa.csv'
B7 = 'fold_B7_aa.csv'
B7A = 'loop_B7_A_aa.csv'
A = 'fold_A_aa.csv'
AC = 'loop_A_C_aa.csv'

totPerResidues = function(listfiles) {
  for (filename in listfiles) {
    d = read.csv(filename)
    d = d[-which(d$pdb =='BEM2_2' | d$pdb=='CLA4Ct' | d$pdb=='NVJ2' | d$pdb=='SKG3N' | d$pdb=='3wyfE00'),]
    if ((exists('tot'))==FALSE) {
      tot = d[,-1]
    } else {
      tot = tot + d[,-1]
    }
  }
  return(tot)
}

totPerFam = function(m) {
  small = apply(m[,which(colnames(m)=='D' | colnames(m)=='N' | colnames(m)=='S' | 
    colnames(m)=='T' | colnames(m)=='C' | colnames(m)=='A' | colnames(m)=='V' | 
    colnames(m)=='P' | colnames(m)=='G')],1,sum)
  polar = apply(m[,which(colnames(m)=='R' | colnames(m)=='K' | colnames(m)=='D' | 
    colnames(m)=='E' | colnames(m)=='Q' | colnames(m)=='N' | colnames(m)=='H' | 
    colnames(m)=='S' | colnames(m)=='T' | colnames(m)=='Y' | colnames(m)=='C' | 
    colnames(m)=='W')],1,sum)
  hydrophobic = apply(m[,which(colnames(m)=='K' | colnames(m)=='H' | colnames(m)=='T' | 
    colnames(m)=='Y' | colnames(m)=='C' | colnames(m)=='W' | colnames(m)=='A' | colnames(m)=='I' | 
    colnames(m)=='L' | colnames(m)=='M' | colnames(m)=='F' | colnames(m)=='V' |colnames(m)=='G')],1,sum)
  tiny = apply(m[,which(colnames(m)=='S' | colnames(m)=='C' | colnames(m)=='A' | colnames(m)=='G')],1,sum)
  aliphatic = apply(m[,which(colnames(m)=='I' | colnames(m)=='L' | colnames(m)=='V')],1,sum)
  aromatic = apply(m[,which(colnames(m)=='H' | colnames(m)=='Y' | colnames(m)=='W' | colnames(m)=='F')],1,sum)
  pos_charged = apply(m[,which(colnames(m)=='R' | colnames(m)=='K' | colnames(m)=='H')],1,sum)
  neg_charged = apply(m[,which(colnames(m)=='D' | colnames(m)=='E')],1,sum)
  return(cbind(small,tiny,polar,hydrophobic,aliphatic,aromatic,pos_charged,neg_charged))
}

error.bar <- function(x, y, upper, lower=upper, length=0.05,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

sem.calculation = function(m) {
  sems = NULL
  for (i in 1:dim(m)[2]){
    sems = c(sems,sd(m[,i])/sqrt(length(m[,i])))
  }
  names(sems) = colnames(m)
  return(sems)
}

mB1B2.aa = totPerResidues(c(B1,B1B2,B2))
mB1B2 = totPerFam(mB1B2.aa)
mB1B2.prot.tot = apply(mB1B2,1,sum)
mB1B2.fam.tot = apply(mB1B2,2,sum)
mB1B2.freq = mB1B2 / mB1B2.prot.tot
mB1B2.tot.freq = mB1B2.fam.tot / sum(mB1B2.fam.tot)
mB1B2.sem = sem.calculation(mB1B2.freq)

mB7.aa = totPerResidues(c(B7))
mB7 = totPerFam(mB7.aa)
mB7.prot.tot = apply(mB7,1,sum)
mB7.fam.tot = apply(mB7,2,sum)
mB7.freq = mB7 / mB7.prot.tot
mB7.tot.freq = mB7.fam.tot / sum(mB7.fam.tot)
mB7.sem = sem.calculation(mB7.freq)

solvent.aa = totPerResidues(c(NB1,B2B3,B3,B3B4,B4,B4B5,B5,B5B6,B6,B6B7,B7A,A,AC))
solvent = totPerFam(solvent.aa)
solvent.prot.tot = apply(solvent,1,sum)
solvent.fam.tot = apply(solvent,2,sum)
solvent.freq = solvent / solvent.prot.tot
solvent.tot.freq = solvent.fam.tot / sum(solvent.fam.tot)
solvent.sem = sem.calculation(solvent.freq)

freq = rbind(mB7.tot.freq,mB1B2.tot.freq,solvent.tot.freq)
sem = rbind(mB7.sem,mB1B2.sem,solvent.sem)
colnames(freq)[which(colnames(freq)=='pos_charged')]='+ charged'
colnames(freq)[which(colnames(freq)=='neg_charged')]='- charged'
colnames(sem) = colnames(freq)

png(filename = "2.30.29.30.4_2_0.inventory.png",width = 1200, height = 800,bg = "white")
par(mai=c(1,1.2,0.35,0.1),mgp=c(4,1,0),cex.axis=1.6,cex.lab=2,family="Bookman Old Style")
barx = barplot(freq,y=c(0,0.35), beside = TRUE,col = c('lightcoral','lightgoldenrod','lightblue'),xaxt = 'n',xlab='',ylab='Frequency')
error.bar(barx,freq,sem)
axis(1, at=seq(2.5,32,4),labels=colnames(freq),cex.axis=1.6,tick=FALSE)
title(xlab='Amino Acid Composition',mgp=c(3.5, 1, 0),cex.lab=2)
legend('topright',bg='gray96',box.col='black',c('β7','β1, β1-β2, β2','solvent side'),col = c('lightcoral','lightgoldenrod','lightblue'),pch=15,cex=1.6,pt.cex=3)
dev.off()