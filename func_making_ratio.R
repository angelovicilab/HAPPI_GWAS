
aa_dat <- data.frame(V1=c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", 
                          "Phe", "Pro", "Pyl", "Ser", "Sec", "Thr", "Trp", "Tyr", "Val", "Asx", "Glx", "Xaa", "Xle"),
                     V2=c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "O", "S", "U", 
                          "T", "W", "Y", "V", "B", "Z", "X", "J"), 
                     V3=c("Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamine", 
                          "Glutamic acid", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", 
                          "Methionine", "Phenylalanine", "Proline", "Pyrrolysine", "Serine", "Selenocysteine", 
                          "Threonine", "Tryptophan", "Tyrosine", "Valine", "Aspartic acid or Asparagine", 
                          "Glutamic acid or Glutamine", "Any amino acid", "Leucine or Isoleucine"),
                     stringsAsFactors = FALSE)

gln_operation <- c("Q/R", "Q/E", "Q/H", "Q/P", "Q/Total", "Q/RQ", "Q/RE", "Q/RH", "Q/RP", "Q/QE", "Q/HP", "Q/QH", "Q/QP", "Q/EH", 
                   "Q/EP", "Q/RQE", "Q/RQH", "Q/RQP", "Q/REH", "Q/REP", "Q/RHP", "Q/QEH", "Q/QEP", "Q/QHP", "Q/EHP", "Q/RQEH", 
                   "Q/RQEP", "Q/RQHP", "Q/REHP", "Q/QEHP", "Q/RQEHP")



maize_operations <- c("IMTXK", "X/IMTXK", "A/X", "X/K", "Z/IMTXK", "Z/I", "Z/K", "I/IMTXK", "K/IMTXK", "M/IMTXK", "M/I", "M/K", 
                      "T/IMTXK", "T/Z", "T/I", "T/K", "ZHPR", "R/ZHPR", "R/Z", "R/P", "X/ZHPR", "Z/ZHPR", "Z/R", "Z/P", "H/ZHPR", 
                      "P/ZHPR", "H/M", "H/Z", "H/X", "IVL", "I/V", "I/IVL", "V/IVL", "Z/V", "LAV", "I/LAV", "A/LAV", "L/V", 
                      "L/LAV","I/L", "L/IVL", "V/LAV", "V/A", "F/FY", "F/Z", "FY", "Y/F", "Y/FY", "G/SG", "S/SG", "A/VA", "E/EL", 
                      "F/FY", "F/FZ", "I/IP", "K/RK", "L/EL", "P/TP", "V/VA", "Z/FZ", "I/LVI", "I/LVIF", "F/LVIF", "XSAVIT", "A/XSAVIT", 
                      "X/XSAVIT", "V/XSAVIT", "T/XSAVIT", "I/XSAVIT", "X/XSAVIT")



generate_gln_ratio <- function(dat){
  if(all(c("Q", "R") %in% colnames(dat))){
    dat$`Q/R` <- dat$Q/dat$R
  }
  if(all(c("Q", "E") %in% colnames(dat))){
    dat$`Q/E` <- dat$Q/dat$E
  }
  if(all(c("Q", "H") %in% colnames(dat))){
    dat$`Q/H` <- dat$Q/dat$H
  }
  if(all(c("Q", "P") %in% colnames(dat))){
    dat$`Q/P` <- dat$Q/dat$P
  }
  if(all(c("Q", "Total") %in% colnames(dat))){
    dat$`Q/Total` <- dat$Q/dat$Total
  }
  if(all(c("Q", "R") %in% colnames(dat))){
    dat$`Q/RQ` <- dat$Q/(dat$R+dat$Q)
  }
  if(all(c("Q", "R", "E") %in% colnames(dat))){
    dat$`Q/RE` <- dat$Q/(dat$R+dat$E)
  }
  if(all(c("Q", "R", "H") %in% colnames(dat))){
    dat$`Q/RH` <- dat$Q/(dat$R+dat$H)
  }
  if(all(c("Q", "R", "P") %in% colnames(dat))){
    dat$`Q/RP` <- dat$Q/(dat$R+dat$P)
  }
  if(all(c("Q", "E") %in% colnames(dat))){
    dat$`Q/QE` <- dat$Q/(dat$Q+dat$E)
  }
  if(all(c("Q", "H", "P") %in% colnames(dat))){
    dat$`Q/HP` <- dat$Q/(dat$H+dat$P)
  }
  if(all(c("Q", "H") %in% colnames(dat))){
    dat$`Q/QH` <- dat$Q/(dat$Q+dat$H)
  }
  if(all(c("Q", "P") %in% colnames(dat))){
    dat$`Q/QP` <- dat$Q/(dat$Q+dat$P)
  }
  if(all(c("Q", "E", "H") %in% colnames(dat))){
    dat$`Q/EH` <- dat$Q/(dat$E+dat$H)
  }
  if(all(c("Q", "E", "P") %in% colnames(dat))){
    dat$`Q/EP` <- dat$Q/(dat$E+dat$P)
  }
  if(all(c("Q", "R", "E") %in% colnames(dat))){
    dat$`Q/RQE` <- dat$Q/(dat$R+dat$Q+dat$E)
  }
  if(all(c("Q", "R", "H") %in% colnames(dat))){
    dat$`Q/RQH` <- dat$Q/(dat$R+dat$Q+dat$H)
  }
  if(all(c("Q", "R", "P") %in% colnames(dat))){
    dat$`Q/RQP` <- dat$Q/(dat$R+dat$Q+dat$P)
  }
  if(all(c("Q", "R", "E", "H") %in% colnames(dat))){
    dat$`Q/REH` <- dat$Q/(dat$R+dat$E+dat$H)
  }
  if(all(c("Q", "R", "E", "P") %in% colnames(dat))){
    dat$`Q/REP` <- dat$Q/(dat$R+dat$E+dat$P)
  }
  if(all(c("Q", "R", "H", "P") %in% colnames(dat))){
    dat$`Q/RHP` <- dat$Q/(dat$R+dat$H+dat$P)
  }
  if(all(c("Q", "E", "H") %in% colnames(dat))){
    dat$`Q/QEH` <- dat$Q/(dat$Q+dat$E+dat$H)
  }
  if(all(c("Q", "E", "P") %in% colnames(dat))){
    dat$`Q/QEP` <- dat$Q/(dat$Q+dat$E+dat$P)
  }
  if(all(c("Q", "H", "P") %in% colnames(dat))){
    dat$`Q/QHP` <- dat$Q/(dat$Q+dat$H+dat$P)
  }
  if(all(c("Q", "E", "H", "P") %in% colnames(dat))){
    dat$`Q/EHP` <- dat$Q/(dat$E+dat$H+dat$P)
  }
  if(all(c("Q", "R", "E", "H") %in% colnames(dat))){
    dat$`Q/RQEH` <- dat$Q/(dat$R+dat$Q+dat$E+dat$H)
  }
  if(all(c("Q", "R", "E", "P") %in% colnames(dat))){
    dat$`Q/RQEP` <- dat$Q/(dat$R+dat$Q+dat$E+dat$P)
  }
  if(all(c("Q", "R", "H", "P") %in% colnames(dat))){
    dat$`Q/RQHP` <- dat$Q/(dat$R+dat$Q+dat$H+dat$P)
  }
  if(all(c("Q", "R", "E", "H", "P") %in% colnames(dat))){
    dat$`Q/REHP` <- dat$Q/(dat$R+dat$E+dat$H+dat$P)
  }
  if(all(c("Q", "E", "H", "P") %in% colnames(dat))){
    dat$`Q/QEHP` <- dat$Q/(dat$Q+dat$E+dat$H+dat$P)
  }
  if(all(c("Q", "R", "E", "H", "P") %in% colnames(dat))){
    dat$`Q/RQEHP` <- dat$Q/(dat$R+dat$Q+dat$E+dat$H+dat$P)
  }
  return(dat)
}



generate_maize_ratio <- function(dat){
  if(all(c("A", "X") %in% colnames(dat))){
    dat$`A/X` <- dat$A/dat$X
  }
  if(all(c("F", "Z") %in% colnames(dat))){
    dat$`F/Z` <- dat$`F`/dat$Z
  }
  if(all(c("H", "M") %in% colnames(dat))){
    dat$`H/M` <- dat$H/dat$M
  }
  if(all(c("H", "X") %in% colnames(dat))){
    dat$`H/X` <- dat$H/dat$X
  }
  if(all(c("H", "Z") %in% colnames(dat))){
    dat$`H/Z` <- dat$H/dat$Z
  }
  if(all(c("I", "L") %in% colnames(dat))){
    dat$`I/L` <- dat$I/dat$L
  }
  if(all(c("I", "V") %in% colnames(dat))){
    dat$`I/V` <- dat$I/dat$V
  }
  if(all(c("L", "V") %in% colnames(dat))){
    dat$`L/V` <- dat$L/dat$V
  }
  if(all(c("M", "I") %in% colnames(dat))){
    dat$`M/I` <- dat$M/dat$I
  }
  if(all(c("M", "K") %in% colnames(dat))){
    dat$`M/K` <- dat$M/dat$K
  }
  if(all(c("P", "R") %in% colnames(dat))){
    dat$`R/P` <- dat$R/dat$P
  }
  if(all(c("Z", "R") %in% colnames(dat))){
    dat$`R/Z` <- dat$R/dat$Z
  }
  if(all(c("T", "I") %in% colnames(dat))){
    dat$`T/I` <- dat$`T`/dat$I
  }
  if(all(c("T", "K") %in% colnames(dat))){
    dat$`T/K` <- dat$`T`/dat$K
  }
  if(all(c("Z", "T") %in% colnames(dat))){
    dat$`T/Z` <- dat$`T`/dat$Z
  }
  if(all(c("V", "A") %in% colnames(dat))){
    dat$`V/A` <- dat$V/dat$A
  }
  if(all(c("K", "X") %in% colnames(dat))){
    dat$`X/K` <- dat$X/dat$K
  }
  if(all(c("Z", "I") %in% colnames(dat))){
    dat$`Z/I` <- dat$Z/dat$I
  }
  if(all(c("Z", "K") %in% colnames(dat))){
    dat$`Z/K` <- dat$Z/dat$K
  }
  if(all(c("Z", "P") %in% colnames(dat))){
    dat$`Z/P` <- dat$Z/dat$P
  }
  if(all(c("Z", "R") %in% colnames(dat))){
    dat$`Z/R` <- dat$Z/dat$R
  }
  if(all(c("Z", "V") %in% colnames(dat))){
    dat$`Z/V` <- dat$Z/dat$V
  }
  
  if(all(c("F", "Y") %in% colnames(dat))){
    dat$`FY` <- dat$`F`+dat$Y
    dat$`Y/F` <- dat$Y/dat$`F`
    dat$`F/FY` <- dat$`F`/dat$`F`+dat$Y
    dat$`Y/FY` <- dat$Y/(dat$`F`+dat$Y)
  }
  
  if(all(c("S", "G") %in% colnames(dat))){
    dat$`G/SG` <- dat$G/(dat$S+dat$G)
    dat$`S/SG` <- dat$S/(dat$S+dat$G)
  }
  
  if(all(c("I", "V", "L") %in% colnames(dat))){
    dat$`IVL` <- dat$I + dat$V + dat$L
    dat$`I/IVL` <- dat$I/dat$`IVL`
    dat$`V/IVL` <- dat$V/dat$`IVL`
    dat$`L/IVL` <- dat$L/dat$`IVL`
  }
  
  if(all(c("L", "A", "V") %in% colnames(dat))){
    dat$`LAV` <- dat$L+dat$A+dat$V
    dat$`L/LAV` <- dat$L/dat$`LAV`
    dat$`A/LAV` <- dat$A/dat$`LAV`
    dat$`V/LAV` <- dat$V/dat$`LAV`
  }
  if(all(c("I", "L", "A", "V") %in% colnames(dat))){
    dat$`I/LAV` <- dat$I/dat$`LAV`
  }
  
  if(all(c("Z", "H", "P", "R") %in% colnames(dat))){
    dat$`ZHPR` <- dat$Z+dat$H+dat$P+dat$R
    dat$`Z/ZHPR` <- dat$Z/dat$`ZHPR`
    dat$`H/ZHPR` <- dat$H/dat$`ZHPR`
    dat$`R/ZHPR` <- dat$R/dat$`ZHPR`
    dat$`P/ZHPR` <- dat$P/dat$`ZHPR`
  }
  if(all(c("X", "Z", "H", "P", "R") %in% colnames(dat))){
    dat$`X/ZHPR` <- dat$X/dat$`ZHPR`
  }
  if(all(c("I", "M", "T", "X", "K") %in% colnames(dat))){
    dat$`IMTXK` <- dat$I+dat$M+dat$`T`+dat$X+dat$K
    dat$`I/IMTXK` <- dat$I/dat$`IMTXK`
    dat$`M/IMTXK` <- dat$M/dat$`IMTXK`
    dat$`T/IMTXK` <- dat$`T`/dat$`IMTXK`
    dat$`X/IMTXK` <- dat$X/dat$`IMTXK`
    dat$`K/IMTXK` <- dat$K/dat$`IMTXK`
  }
  if(all(c("Z", "I", "M", "T", "X", "K") %in% colnames(dat))){
    dat$`Z/IMTXK` <- dat$Z/dat$`IMTXK`
  }
  if(all(c("A", "V") %in% colnames(dat))){
    dat$`A/VA` <- dat$A/(dat$A+dat$V)
    dat$`V/VA` <- dat$V/(dat$A+dat$V)
  }
  if(all(c("E", "L") %in% colnames(dat))){
    dat$`E/EL` <- dat$E/(dat$E+dat$L)
  }
  if(all(c("F", "Y") %in% colnames(dat))){
    dat$`F/FY` <- dat$`F`/(dat$`F`+dat$Y)
  }
  if(all(c("F", "Z") %in% colnames(dat))){
    dat$`F/FZ` <- dat$`F`/(dat$`F`+dat$Z)
  }
  if(all(c("F", "Y") %in% colnames(dat))){
    dat$`I/IP` <- dat$I/(dat$I+dat$P)
  }
  if(all(c("F", "Y") %in% colnames(dat))){
    dat$`K/RK` <- dat$I/(dat$I+dat$P)
  }
  if(all(c("E", "L") %in% colnames(dat))){
    dat$`L/EL` <- dat$L/(dat$E+dat$L)
  }
  if(all(c("E", "L") %in% colnames(dat))){
    dat$`L/EL` <- dat$L/(dat$E+dat$L)
  }
  if(all(c("P", "T") %in% colnames(dat))){
    dat$`P/TP` <- dat$P/(dat$P+dat$`T`)
  }
  if(all(c("Z", "F") %in% colnames(dat))){
    dat$`Z/FZ` <- dat$Z/(dat$`F`+dat$Z)
  }
  if(all(c("I", "L", "V") %in% colnames(dat))){
    dat$`I/LVI` <- dat$I/(dat$L+dat$V+dat$I)
  }
  if(all(c("I", "L", "V", "F") %in% colnames(dat))){
    dat$`I/LVIF` <- dat$I/(dat$L+dat$V+dat$I+dat$`F`)
    dat$`F/LVIF` <- dat$`F`/(dat$L+dat$V+dat$I+dat$`F`)
  }
  if(all(c("X", "S", "A", "V", "I", "T") %in% colnames(dat))){
    dat$`XSAVIT` <- dat$X+dat$S+dat$A+dat$V+dat$I+dat$`T`
    dat$`X/XSAVIT` <- dat$X/(dat$X+dat$S+dat$A+dat$V+dat$I+dat$`T`)
    dat$`S/XSAVIT` <- dat$S/(dat$X+dat$S+dat$A+dat$V+dat$I+dat$`T`)
    dat$`A/XSAVIT` <- dat$A/(dat$X+dat$S+dat$A+dat$V+dat$I+dat$`T`)
    dat$`V/XSAVIT` <- dat$V/(dat$X+dat$S+dat$A+dat$V+dat$I+dat$`T`)
    dat$`I/XSAVIT` <- dat$I/(dat$X+dat$S+dat$A+dat$V+dat$I+dat$`T`)
    dat$`T/XSAVIT` <- dat$`T`/(dat$X+dat$S+dat$A+dat$V+dat$I+dat$`T`)
  }
  
  return(dat)
}
