library(seqinr)
library(ape)
GoIATPAM<-function(contig,insertlocation,closedistance=50,primersnumber=5){
  contig[is.na(contig)]<-"N"
  ATPAM<-gregexpr(pattern ='a.{5,7}t.{4,6}.gg',paste(contig,collapse=""))[[1]]

  ATPAM[tail(which(insertlocation-ATPAM>closedistance),primersnumber)]
}


CRISPRprimercall<-function(contig,GoIATPAM,gRNA,tRNA,gRNApLength=23,tRNApLength=21){
  primers<-list()
  insertprecise<-list()
  stickyend<-list()
  for(p in c(1:length(GoIATPAM))){
    seq<-contig[GoIATPAM[p]:(GoIATPAM[p]+20)]
    for(s in 7:9){
      if(seq[s]=="t"){
        seq[s]="u"
      }
    }
    for(s in 7:9){
      if(sum(seq=="u")>1&seq[s]=="u"){
        seq[s]<-"t"
      }
    }
    for(s in (which(seq=="u")+6):(which(seq=="u")+8)){
      if(seq[s]=="g"&seq[s+1]=="g"){
        seq[s]="x"
        seq[s+1]="x"
        break
      }
    }
    seq<-seq[1:which(seq=="x")[2]]
    PF<-paste0(paste0(seq[1:(length(seq)-3)],collapse=""),substring(gRNA,1,gRNApLength))
    bp20<-contig[(GoIATPAM[p]-20+length(seq)-3):(GoIATPAM[p]-1+length(seq)-3)]
    compseq<-rev(comp(bp20))[(length(seq)-3-which(seq=="u")+1):20]
    compseq[which(seq=="u")]<-"u"
    PR<-paste0(paste0(compseq,collapse=""),substring(tRNA,1,tRNApLength))
    primers[[p]]<-c(PF,PR)
    insertprecise[[p]]<-c(GoIATPAM[p]-20+length(seq)-3,GoIATPAM[p]-1+length(seq)-3)
    stickyend[[p]]<-c(paste(seq[1:which(seq=="u")],collapse=""),paste(compseq[1:which(seq=="u")],collapse=""))
  }
  return(list(primers=primers,insertsite=insertprecise,stickyend=stickyend))
}


genome<-read.fasta("data/genome/GCA_021066465.1_ASM2106646v1_genomic.fna")
gff<-read.gff("data/genome/GCA_021066465.1_ASM2106646v1_genomic.gff")
GoI<-"LI328DRAFT_134294"
genelocation<-gff[gff$type=="gene"&grepl(GoI,gff$attributes),]

insertlocation<-c(as.character(genelocation$seqid),ifelse(genelocation$strand=="+",genelocation$start,genelocation$end))


##PU3-tRNA-PS1-gRNA-tRNA-TU3

tRNA<-"GCATCATTGGTCTAGTGGTAGAATTCATCGTTGCCATCGATGAGGCCCGTGTTCGATTCACGGATGATGCA"
gRNA<-"GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT"

tRNAprimerfrag<-seqinr::s2g(seqinr::s2c(tRNA)[(nchar(tRNA)-21):nchar(tRNA)])
gRNAprimerfrag<-seqinr::s2g(seqinr::s2c(gRNA)[(nchar(gRNA)-23):nchar(gRNA)])

#Find  upstreamAT-PAM (NGG)
contig<-genome[[insertlocation[1]]]
upstreamGoIATPAM<-GoIATPAM(contig,as.numeric(insertlocation[2]))
upstreamprimers<-CRISPRprimercall(contig,upstreamGoIATPAM,gRNA,tRNA,gRNApLength=23,tRNApLength=21)

#Find  downstreamAT-PAM (NGG)
contig<-rev(comp(genome[[insertlocation[1]]]))
downstreamGoIATPAM<-GoIATPAM(contig,length(contig)-as.numeric(insertlocation[2]))
downstreamprimers<-CRISPRprimercall(contig,downstreamGoIATPAM,gRNA,tRNA,gRNApLength=23,tRNApLength=21)

##pick one

#
