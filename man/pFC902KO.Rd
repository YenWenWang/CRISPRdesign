% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pFC902.R
\name{pFC902KO}
\alias{pFC902KO}
\title{pFC902 KO design}
\usage{
pFC902KO(
  contig,
  siteofinterest,
  offset = 0,
  upstream = T,
  primersnumber = 5,
  gRNApLength = 23,
  tRNApLength = 21
)
}
\arguments{
\item{contig}{Sequence of contig that include site of interest}

\item{siteofinterest}{Approximate site for double strand break}

\item{offset}{give some space between the two double strand break sites}

\item{upstream}{Find site upstream or downstream}

\item{primersnumber}{How many primers you want}
}
\value{
a list of suggested primers, additional primers and primers required for pFC330~pFC333 cloning
}
\description{
design primers for knockout using pFC902
}
\details{
TBD
}
\examples{
genome<-seqinr::read.fasta("data/genome/GCA_021066465.1_ASM2106646v1_genomic.fna")
gff<-ape::read.gff("data/genome/GCA_021066465.1_ASM2106646v1_genomic.gff")
GoI<-"LI328DRAFT_134294"
genelocation<-gff[gff$type=="gene"&grepl(GoI,gff$attributes),]
contig<-genome[[as.character(genelocation$seqid)]]
siteofinterest<-ifelse(genelocation$strand=="+",genelocation$start,genelocation$end)
pFC902KO(contig,siteofinterest,50)
}
\references{
Nødvig, et al. 2018. Efficient oligo nucleotide mediated CRISPR-Cas9 gene editing in Aspergilli. Fungal Genet. Biol. 115:78-89.
}
