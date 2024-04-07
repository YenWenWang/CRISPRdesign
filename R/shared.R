#' Common insert data
#'
#' @description
#' Sequence for some common inserts
#'
#' @details
#' include pTCU1 from Trichoderma reesei
promoter_pTCU1<-"CTGTGTGGCATCACTCATGTTCTGGATGTGCCGAGCACGCACTATAGAGTGAGCGCCAGCCATGGCGTTGCGCATTTGGTCTCTTGATAGAAAAGGCAGGACGGGTGACTCTGTGCCAGATGGTTGCGGCGCCGAGTGTTGGTGATTGGCAGCGTTATGTTGCAAAGGAGCTTTATCGGGCAGCACAGAGAACTTAATATAAGCATGACGAGCAGCCAGATAAGTTCAATACCCAAGAAATTCTAGGGCTCTAGTTATGTTGTCCACTTGGGAGGCTTGGAGGGATGGTCCCCCCCCCCCCCCTTAAAGTAATTCGCGGCTGTTGTATCGGTGCGCTTTATCCGCCAACGGTTTTGATAGTTGCACTTGTATCTGCTGGAATGCTTCCACAATGCTGTAGTCAGTGATACAAGGTTAATGATGTCGCTTGATGCGGATTCTTGAAGGTGGGGAACTACTGCAGGTGAGGACATATGAGCAAGTTTGCGGGAGATTATGAAGATCGCAGGTTGGCTGACTTGACCACCCTTATCAGGCACATTTGATTCGGCCTATATATGTGTGAGATTTATCTCAGCAATGCTCAGACAACTCCTCTGAATGGTTTGGCGATTTTGCTTTGGCTTGATCATCGCTTGTCTTTGCAGTCACCCTACTTCACAGGACTCCCAAGGAGTGCATCAATCCACAAGAGCCTACTGCCAAATCCAGGAGTTCAATCTATACGACCTGGTTGATACGACA"


#' PAM site for USER
#'
#' @description Find PAM sites with A and T upstream for Uracil Specific Excision Reaction
#'
#' @param contig
#' Sequence of contig that include site of interest
#'
#' @param siteofinterest
#' Approximate site for double strand break
#'
#' @param upstream
#' Find site upstream or downstream
#'
#' @param primersnumber
#' How many PAMs you want
#'
#' @return exact
#' @details
#' an internal function TBD
#'
#' @examples
#' USERPAM(contig,site-50)
#'
USERPAM<-function(contig,siteofinterest,upstream=T,primersnumber=5){
  contig[is.na(contig)]<-"N"
  USERPAM<-gregexpr(pattern ='[Aa].{5,7}[Tt].{4,6}.gg',paste(contig,collapse=""))[[1]]
  if(upstream){
    USERPAM[rev(tail(which(siteofinterest-USERPAM>0),primersnumber))]
  }else{
    USERPAM[head(which(USERPAM-siteofinterest>0),primersnumber)]
  }
}


#' Call CRISPR primers
#'
#' @description The function finds the PAM sites and design primers to incorporate cleavage sites.
#'
#' @param contig
#' Sequence of contig that include site of interest
#'
#' @param siteofinterest
#' Approximate site for double strand break
#'
#' @param PFrecogsite
#' Sequence attached to forward strand of the target recognition seq.
#'
#' @param PRrecogsite
#' Sequence attached to reverse strand of the target recognition seq.
#'
#' @param upstream
#' Find site upstream or downstream
#'
#' @param primersnumber
#' How many primer pairs you want
#'
#' @return a list of primers, locations and sticky ends
#' @details
#' TBD
#'
#' @export
#'
#' @examples
#' CRISPRUSERprimercall(contig,site-50,PFrecogsite,PRrecogsite)
#'
CRISPRUSERprimercall<-function(contig,siteofinterest,PFrecogsite,PRrecogsite,upstream=T,primersnumber=5){
  USERPAMsites<-USERPAM(contig,siteofinterest,upstream,primersnumber)
  primers<-list()
  siteprecise<-list()
  stickyend<-list()
  for(p in c(1:primersnumber)){
    seq<-contig[USERPAMsites[p]:(USERPAMsites[p]+20)]
    for(s in 7:9){
      if(seq[s]%in%c("T","t")){
        seq[s]="u"
      }
    }
    for(s in 7:9){
      if(sum(seq=="u")>1&seq[s]=="u"){
        seq[s]<-"t"
      }
    }
    for(s in (which(seq=="u")+6):(which(seq=="u")+8)){
      if(seq[s]%in%c("G","g")&seq[s+1]%in%c("G","g")){
        seq[s]="x"
        seq[s+1]="x"
        break
      }
    }
    seq<-seq[1:which(seq=="x")[2]]
    PF<-paste0(paste0(seq[1:(length(seq)-3)],collapse=""),PFrecogsite)
    bp20<-contig[(USERPAMsites[p]-20+length(seq)-3):(USERPAMsites[p]-1+length(seq)-3)]
    compseq<-rev(seqinr::comp(bp20))[(length(seq)-3-which(seq=="u")+1):20]
    compseq[which(seq=="u")]<-"u"
    PR<-paste0(paste0(compseq,collapse=""),PRrecogsite)
    primers[[p]]<-c(PF,PR)
    siteprecise[[p]]<-c(USERPAMsites[p]-20+length(seq)-3,USERPAMsites[p]-1+length(seq)-3)
    stickyend[[p]]<-c(paste(seq[1:which(seq=="u")],collapse=""),paste(compseq[1:which(seq=="u")],collapse=""))
  }
  return(list(primers=primers,site=siteprecise,stickyend=stickyend))
}


#' USER site at ends
#'
#' @description Find sites Uracil Specific Excision Reaction at 5' and 3' ends of a sequence
#'
#' @param sequence
#' Sequence of of interest
#'
#' @return the start and end of the two USER sites
#' @details
#' TBD
#'
#' @examples
#' USERends(pTCU1)
#'
USERends<-function(sequence){
  first<-regexpr(pattern ='[Aa].{5,7}[Tt]',sequence)
  last<-regexpr(pattern ='[Aa].{5,7}[Tt]',paste0(rev(seqinr::comp(seqinr::s2c(sequence))),collapse = ""))
  return(list(first=c(first,first+attr(first,"match.length")-1),
              last=c(nchar(pTCU1)-last-attr(last,"match.length")+2,nchar(pTCU1)-last+1)))
}
