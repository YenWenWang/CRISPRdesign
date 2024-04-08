#' pFC902 KO design
#'
#' @description design primers for knockout using pFC902
#'
#' @param contig
#' Sequence of contig that include site of interest
#'
#' @param siteofinterest
#' Approximate site for double strand break
#'
#' @param offset
#' give some space between the two double strand break sites
#'
#' @param upstream
#' Find site upstream or downstream
#'
#' @param primersnumber
#' How many primers you want
#'
#' @return a list of suggested primers, additional primers and primers required for pFC330~pFC333 cloning
#' @details
#' TBD
#'
#' @references
#' Nødvig, et al. 2018. Efficient oligo nucleotide mediated CRISPR-Cas9 gene editing in Aspergilli. Fungal Genet. Biol. 115:78-89.
#'
#' @examples
#' genome<-seqinr::read.fasta("data/genome/GCA_021066465.1_ASM2106646v1_genomic.fna")
#' gff<-ape::read.gff("data/genome/GCA_021066465.1_ASM2106646v1_genomic.gff")
#' GoI<-"LI328DRAFT_134294"
#' genelocation<-gff[gff$type=="gene"&grepl(GoI,gff$attributes),]
#' contig<-genome[[as.character(genelocation$seqid)]]
#' siteofinterest<-ifelse(genelocation$strand=="+",genelocation$start,genelocation$end)
#' pFC902KO(contig,siteofinterest,50)

pFC902KO<-function(contig,siteofinterest,offset=0,upstream=T,primersnumber=5,gRNApLength=23,tRNApLength=21){
  PFrecogsite<-substring(pFC902_gRNA,1,gRNApLength)
  PRrecogsite<-toupper(substring(paste0(rev(seqinr::comp(seqinr::s2c(pFC902_tRNA))),collapse = ""),1,tRNApLength))

  site1<-CRISPRUSERprimercall(contig,siteofinterest,PFrecogsite,PRrecogsite,upstream=T,primersnumber=5)

  compcontig<-rev(seqinr::comp(contig))
  site2<-CRISPRUSERprimercall(compcontig,length(contig)-siteofinterest,PFrecogsite,PRrecogsite,upstream=T,primersnumber=5)
  site2$site<-lapply(site2$site,function(x)length(contig)-x+1)

  if(any(duplicated(c(unlist(site1$stickyend),unlist(site2$stickyend),
               gsub("u.*","u",tolower(c(primer_CSN438,primer_CSN790))))))){
    return(list(recommended="some sticky ends are duplicated; please hand pick",
                site1=site1,site2=site2,required=c(primer_CSN438,primer_CSN790)))
  }

  pair1<-c(primer_CSN438,site1$primers[[1]][2])
  pair2<-c(site1$primers[[1]][1],site2$primers[[1]][2])
  pair3<-c(site2$primers[[1]][1],primer_CSN790)

  warning("please check if the primers are too similar, Tm or have hairpin etc.")
  return(list(recommended=list(pair1,pair2,pair3),
              site1=site1,site2=site2,required=c(primer_CSN438,primer_CSN790)))

}




#' pFC902 insert design
#'
#' @description design primers to introduce an insert into genome using pFC902
#'
#' @param contig
#' Sequence of contig that include site of interest
#'
#' @param siteofinterest
#' Insert site
#'
#' @param insert
#' Sequence of insert
#'
#' @param direction
#' + or - (check if your insert is the same direction as your gene of interest for example)
#'
#' @param HMR1primers
#' A vector contains two sequences that contains the whole HMR1 region \
#' (I think it's easier to use NCBI primer finder for this to avoid specificity issues etc.)
#'
#' @param HMR2primers
#' A vector contains two sequences that contains the whole HMR2 region \
#' (I think it's easier to use NCBI primer finder for this to avoid specificity issues etc.)
#'
#' @param offset
#' Give some space between the two double strand break sites
#'
#' @param upstream
#' Find site upstream or downstream
#'
#' @param primersnumber
#' How many primers you want for double strand break
#'
#' @return a list of suggested primers, additional primers and primers required for pFC330~pFC333 cloning
#' @details
#' TBD
#'
#' @references
#' Nødvig, et al. 2018. Efficient oligo nucleotide mediated CRISPR-Cas9 gene editing in Aspergilli. Fungal Genet. Biol. 115:78-89.
#'
#' @examples
#' genome<-seqinr::read.fasta("data/genome/GCA_021066465.1_ASM2106646v1_genomic.fna")
#' gff<-ape::read.gff("data/genome/GCA_021066465.1_ASM2106646v1_genomic.gff")
#' GoI<-"LI328DRAFT_134294"
#' genelocation<-gff[gff$type=="gene"&grepl(GoI,gff$attributes),]
#' contig<-genome[[as.character(genelocation$seqid)]]
#' siteofinterest<-ifelse(genelocation$strand=="+",genelocation$start,genelocation$end)
#' pFC902insert(contig,siteofinterest,pTCU1,"+",50) TBD

pFC902insert<-function(contig,siteofinterest,insert,direction,HMR1primers,HMR2primers,
                       offset=50,upstream=T,HMRlength=1500,
                       primersnumber=5,insertpLength=21,Paf_U3pLength=23,gRNApLength=23,tRNApLength=21){
  PFrecogsite<-substring(pFC902_gRNA,1,gRNApLength)
  PRrecogsite<-toupper(substring(paste0(rev(seqinr::comp(seqinr::s2c(pFC902_tRNA))),collapse = ""),1,tRNApLength))



  site1<-CRISPRUSERprimercall(contig,siteofinterest,PFrecogsite,PRrecogsite,upstream=T,primersnumber=5)

  compcontig<-rev(seqinr::comp(contig))
  site2<-CRISPRUSERprimercall(compcontig,length(contig)-siteofinterest,PFrecogsite,PRrecogsite,upstream=T,primersnumber=5)
  site2$site<-lapply(site2$site,function(x)length(contig)-x+1)

  ## FIND insert USER sites
  insertUSERsites<-USERends(insert)
  insertprimerF<-seqinr::s2c(tolower(substring(insert,insertUSERsites$first[1],insertUSERsites$first[1]+insertpLength-1)))
  insertprimerF[insertUSERsites$first[2]-insertUSERsites$first[1]+1]<-"u"
  insertprimerF<-paste0(insertprimerF,collapse = "")
  insertprimerR<-rev(seqinr::comp(seqinr::s2c(substring(insert,insertUSERsites$last[2]-insertpLength+1,insertUSERsites$last[2]))))
  insertprimerR[insertUSERsites$last[2]-insertUSERsites$last[1]+1]<-"u"
  insertprimerR<-paste0(insertprimerR,collapse = "")

  ## FIND Paf_U3 USER sites
  Paf_U3USERsites<-USERends(promoter_Paf_U3)

  ## glue USER sites onto the HMR1 primers.
  HMR1primerF<-paste0(sub("U.*","U",primer_CSN438),HMR1primers[1]) # this will allow the ligation of pFC330~333

  stickyendtemp<-seqinr::comp(rev(seqinr::s2c(tolower(substring(insert,1,insertUSERsites$first[2])))))
  stickyendtemp[length(stickyendtemp)-insertUSERsites$first[1]+1]<-"u"
  HMR1primerR<-paste0(paste0(stickyendtemp,collapse = ""),HMR1primers[2]) # this ligates to insert

  ## glue USER sites onto the HMR2 primers.
  HMR2primerF<-seqinr::s2c(paste0(tolower(substring(insert,insertUSERsites$last[1],nchar(insert))),HMR2primers[1])) # this ligates to insert
  HMR2primerF[insertUSERsites$last[2]-insertUSERsites$last[1]+1]<-"u"
  HMR2primerF<-paste0(HMR2primerF,collapse="")

  stickyendtemp<-seqinr::comp(rev(seqinr::s2c(tolower(substring(promoter_Paf_U3,1,Paf_U3USERsites$first[2])))))
  stickyendtemp[length(stickyendtemp)-Paf_U3USERsites$first[1]+1]<-"u"
  HMR2primerR<-paste0(paste0(stickyendtemp,collapse = ""),HMR2primers[2]) # this ligates to insert paf U3

  output<-list(site1=site1,site2=site2,required=c(primer_CSN790))

  if(any(duplicated(c(unlist(site1$stickyend),unlist(site2$stickyend),
                      gsub("u.*","u",tolower(c(insertprimerF,insertprimerR,HMR1primerF,HMR1primerR,HMR2primerF,HMR2primerR))),
                      gsub("u.*","u",tolower(c(primer_CSN790))))))){
    output$recommended="some sticky ends are duplicated; please hand pick"
    return(output)
  }

  pair1F<-seqinr::s2c(substring(promoter_Paf_U3,Paf_U3USERsites$first[1],Paf_U3USERsites$first[1]+Paf_U3pLength-1))
  pair1F[Paf_U3USERsites$first[2]-Paf_U3USERsites$first[1]+1]<-"U"

  pair1<-c(paste0(pair1F,collapse=""),site1$primers[[1]][2])
  pair2<-c(site1$primers[[1]][1],site2$primers[[1]][2])
  pair3<-c(site2$primers[[1]][1],primer_CSN790)

  warning("please check if the primers are too similar, Tm or have hairpin etc.")
  output$recommended<-list(pair1,pair2,pair3)
  output$insertprimers<-c(insertprimerF,insertprimerR)
  output$HMRprimers<-list(HMR1=c(HMR1primerF,HMR1primerR),HMR2=c(HMR2primerF,HMR2primerR))
  output$protocol<-c(paste("PCR insert DNA with primers",paste(output$insertprimers,collapse = " ")),
                     paste("PCR target DNA with primers",paste(output$HMRprimers$HMR1,collapse = " ")),
                     paste("PCR target DNA with primers",paste(output$HMRprimers$HMR2,collapse = " ")),
                     paste("PCR pFC902 with primers",paste(output$recommended[[1]],collapse = " ")),
                     paste("PCR pFC902 with primers",paste(output$recommended[[2]],collapse = " ")),
                     paste("PCR pFC902 with primers",paste(output$recommended[[3]],collapse = " ")),
                     "digest pFC330~333 with Nt. BbvCI.",
                     "ligate everything together")
  return(output)

}


#' pFC902 tRNA
#'
#' @description
#' pFC902 tRNA
pFC902_tRNA<-"GCATCATTGGTCTAGTGGTAGAATTCATCGTTGCCATCGATGAGGCCCGTGTTCGATTCACGGATGATGCA" #forward

#' pFC902 gRNA
#'
#' @description
#' pFC902 gRNA
pFC902_gRNA<-"GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT" #forward

#' primer CSN438
#'
#' @description
#' primer xc
primer_CSN438<-"GGGTTTAAUGATCACATAGATGCTCGGTTGACA"

#' primer CSN790
#'
#' @description
#' primer CSN790
primer_CSN790<-"GGTCTTAAUACCCTGAGAAGATAGATGTGAATGTG"

#' promoter Paf_U3
#'
#' @description
#' promoter Paf_U3
promoter_Paf_U3<-"GATCACATAGATGCTCGGTTGACAGGACCACCTATAGCAAGTACTTTGTAGGAGGTTCTCCTAACGCTTGGCTTCAAATCCAATCAGGATATCAGTGAGTTACTGTTTGTCATCGCATCCTCATGAAGTCGCTCAACAGCCACTGAGAAGCAAATATTTGGGAGACCATCCCTGATGTTGAAATTTTACGCTGGAGCCCATTCACCGGTGAGCTCGGGGACAGTCTGGCCAGTGGAGCGGAAAATCTTCAACTAACTCTGATTGGTACCACAGGTAGCCCCACAGAACAAACCGAGCAAACATGAAAATTTTCGCTTGAGGTTAGCGCACTCGCTAGCGCCTGCCCGCAAATATAAAAGGCCCCGAAATCCCATGCATTTTGGAAATCAAGCGACTCTACGTATGTCATCCCCCTCCCTACCTTTGGGAAAACCTCTCAACTGCAAGTCAGAACATTTTGCTAACAGC"

#' terminator Taf_U3
#'
#' @description
#' terminator Taf_U3
terminator_Taf_U3<-"TGTCAATTTTTACACTTGATGTTGGTGTAATCAGGCGATCAGAGTTTGCTCAAATTCTCTCTTGTTTTAGATTCTGAAGAATATAATTTTTACATGACTAGTTCATTTTACTTCCAAGTAGACTATTTGTAGATTGAGATTGCCGTCAACTTTATTCATAGAATCACATTAAATTGCTCATTGACCTGCATGTCTATGCACACATTCACATCTATCTTCTCAGGGT"
