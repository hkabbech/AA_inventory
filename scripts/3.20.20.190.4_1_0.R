setwd(paste(getwd(),'/superfamilies/3.20.20.190.4_1_0/tables_old/',sep=''))
B1A1 = 'loop_B1_A1_aa.csv'
A1B2 = 'loop_A1_B2_aa.csv'
B2A2 = 'loop_B2_A2_aa.csv'
A2B3 = 'loop_A2_B3_aa.csv'
B3A3 = 'loop_B3_A3_aa.csv'
A3B4 = 'loop_A3_B4_aa.csv'

totPerResidues = function(listfiles) {
  for (filename in listfiles) {
    d = read.csv(filename)
    d = d[-which(d$pdb=='A0A0D1CR16_USTMA' | d$pdb=='A0CKY8_PARTE' | d$pdb=='A3LSB2_PICST' | d$pdb=='A4S165_OSTLU' | 
                 d$pdb=='A5DAA1_PICGU' | d$pdb=='A5E691_LODEL' | d$pdb=='A6RBI2_AJECN' | d$pdb=='A7T501_NEMVE' | 
                 d$pdb=='A7TR21_VANPO' | d$pdb=='A7UXE0_NEUCR' | d$pdb=='A9RJG9_PHYPA' | d$pdb=='A9SXU1_PHYPA' | 
                 d$pdb=='B6C9K2_SOLLC' | d$pdb=='E9QD74_DANRE' | d$pdb=='G4UMR8_NEUT9' | d$pdb=='H2SZT5_TAKRU' | 
                 d$pdb=='H2USZ2_TAKRU' | d$pdb=='H3CYP8_TETNG' | d$pdb=='H3D277_TETNG' | d$pdb=='H3D2C1_TETNG' | 
                 d$pdb=='H3DK53_TETNG' | d$pdb=='M1BNL7_SOLTU' | d$pdb=='M1CL37_SOLTU' | d$pdb=='M4FG24_BRARP' | 
                 d$pdb=='PLC1_YEAST' | d$pdb=='PLCD1_ARATH' | d$pdb=='PLCD3_HUMAN' | d$pdb=='PLCD4_ARATH' | 
                 d$pdb=='PLCD5_ARATH' | d$pdb=='PLCD8_ARATH' | d$pdb=='PLCD9_ARATH' | d$pdb=='PLCG2_MOUSE' | 
                 d$pdb=='PLCL1_HUMAN' | d$pdb=='PLCL2_MOUSE' | d$pdb=='PLC_DICDI' | d$pdb=='Q0UVV5_PHANO' | 
                 d$pdb=='Q21754_CAEEL' | d$pdb=='Q236T7_TETTS' | d$pdb=='Q2GUM5_CHAGB' | d$pdb=='Q2HB23_CHAGB' | 
                 d$pdb=='Q2UPV6_ASPOR' | d$pdb=='Q4N7P5_THEPA' | d$pdb=='Q4WX21_ASPFU' | d$pdb=='Q5BFL6_EMENI' | 
                 d$pdb=='Q5KJZ5_CRYNJ' | d$pdb=='Q6BZ51_DEBHA' | d$pdb=='Q6CD33_YARLI' | d$pdb=='Q6FY00_CANGA' | 
                 d$pdb=='Q75C92_ASHGO' | d$pdb=='Q7SE08_NEUCR' | d$pdb=='Q8IA75_CAEEL' | d$pdb=='Q75IL8_ORYSJ' |
                 d$pdb=='PLCD6_ARATH'),]
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
                           colnames(m)=='T' | colnames(m)=='C' | colnames(m)=='A' | 
                           colnames(m)=='V' | colnames(m)=='P' | colnames(m)=='G')],1,sum)
  polar = apply(m[,which(colnames(m)=='R' | colnames(m)=='K' | colnames(m)=='D' | 
                           colnames(m)=='E' | colnames(m)=='Q' | colnames(m)=='N' | 
                           colnames(m)=='H' | colnames(m)=='S' | colnames(m)=='T' |
                           colnames(m)=='Y' | colnames(m)=='C' | colnames(m)=='W')],1,sum)
  hydrophobic = apply(m[,which(colnames(m)=='K' | colnames(m)=='H' | colnames(m)=='T' | 
                                 colnames(m)=='Y' | colnames(m)=='C' | colnames(m)=='W' | 
                                 colnames(m)=='A' | colnames(m)=='I' | colnames(m)=='L' |
                                 colnames(m)=='M' | colnames(m)=='F' | colnames(m)=='V' |
                                 colnames(m)=='G')],1,sum)
  tiny = apply(m[,which(colnames(m)=='S' | colnames(m)=='C' | colnames(m)=='A' | colnames(m)=='G')],1,sum)
  aliphatic = apply(m[,which(colnames(m)=='I' | colnames(m)=='L' | colnames(m)=='V')],1,sum)
  aromatic = apply(m[,which(colnames(m)=='H' | colnames(m)=='Y' | colnames(m)=='W' | colnames(m)=='F')],1,sum)
  pos_charged = apply(m[,which(colnames(m)=='R' | colnames(m)=='K' | colnames(m)=='H')],1,sum)
  neg_charged = apply(m[,which(colnames(m)=='D' | colnames(m)=='E')],1,sum)
  fam = cbind(small,tiny,polar,hydrophobic,aliphatic,aromatic,pos_charged,neg_charged)
  return(fam)
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

mB1A1.aa = totPerResidues(c(B1A1))
mB1A1 = totPerFam(mB1A1.aa)
mB1A1.prot.tot = apply(mB1A1,1,sum)
mB1A1.fam.tot = apply(mB1A1,2,sum)
mB1A1.freq = mB1A1 / mB1A1.prot.tot
mB1A1.tot.freq = mB1A1.fam.tot / sum(mB1A1.fam.tot)
mB1A1.sem = sem.calculation(mB1A1.freq)

mB2A2.aa = totPerResidues(c(B2A2))
mB2A2 = totPerFam(mB2A2.aa)
mB2A2.prot.tot = apply(mB2A2,1,sum)
mB2A2.fam.tot = apply(mB2A2,2,sum)
mB2A2.freq = mB2A2 / mB2A2.prot.tot
mB2A2.tot.freq = mB2A2.fam.tot / sum(mB2A2.fam.tot)
mB2A2.sem = sem.calculation(mB2A2.freq)

mB3A3.aa = totPerResidues(c(B3A3))
mB3A3 = totPerFam(mB3A3.aa)
mB3A3.prot.tot = apply(mB3A3,1,sum)
mB3A3.fam.tot = apply(mB3A3,2,sum)
mB3A3.freq = mB3A3 / mB3A3.prot.tot
mB3A3.tot.freq = mB3A3.fam.tot / sum(mB3A3.fam.tot)
mB3A3.sem = sem.calculation(mB3A3.freq)

solvent.aa = totPerResidues(c(A1B2,A2B3,A3B4))
solvent = totPerFam(solvent.aa)
solvent.prot.tot = apply(solvent,1,sum)
solvent.fam.tot = apply(solvent,2,sum)
solvent.freq = solvent / solvent.prot.tot
solvent.tot.freq = solvent.fam.tot / sum(solvent.fam.tot)
solvent.sem = sem.calculation(solvent.freq)

freq = rbind(mB1A1.tot.freq,mB2A2.tot.freq,mB3A3.tot.freq,solvent.tot.freq)
sem = rbind(mB1A1.sem,mB2A2.sem,mB3A3.sem,solvent.sem)
colnames(freq)[which(colnames(freq)=='pos_charged')]='+ charged'
colnames(freq)[which(colnames(freq)=='neg_charged')]='- charged'
colnames(sem) = colnames(freq)

png(filename = "3.20.20.190.4_1_0.inventory.png",width = 1200, height = 800,bg = "white")
par(mai=c(1,1.2,0.35,0.1),mgp=c(4,1,0),cex.axis=1.6,cex.lab=2,family="Bookman Old Style")
barx = barplot(freq,y=c(0,0.32), beside = TRUE,col = c('lightcoral','lightgoldenrod','lightgreen','lightblue'),xaxt = 'n',xlab='',ylab='Frequency')
error.bar(barx,freq,sem)
axis(1, at=seq(3,40,5),labels=colnames(freq),cex.axis=1.6,tick=FALSE)
title(xlab='Amino Acid Composition',mgp=c(3.5, 1, 0),cex.lab=2)
legend('topright',bg='gray96',box.col='black',c('β1-α1','β2-α2','β3-α3','solvent side'),col = c('lightcoral','lightgoldenrod','lightgreen','lightblue'),pch=15,cex=1.6,pt.cex=3)
dev.off()