pvaluelegend <- function(stat,pval){
  rp2 = vector('expression',2)
          
  rp2[1] = substitute(expression("Pearson's r"== MYVALUE),
              list(MYVALUE = format(stat, dig=2)))[2]
  rp2[2] = substitute(expression(italic(P)== MYOTHERVALUE),
              list(MYOTHERVALUE = (ifelse(pval>0.0001,format(pval, digits = 3), "0.0001"))))[2]

  legend("topleft", legend = rp2, bty = 'n', text.col="black", cex=1.3) 
}

#pvaluelegend(cortest1$estimate,cortest1$p.value)

