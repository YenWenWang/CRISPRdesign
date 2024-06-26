% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shared.R
\name{CRISPRUSERprimercall}
\alias{CRISPRUSERprimercall}
\title{Call CRISPR primers}
\usage{
CRISPRUSERprimercall(
  contig,
  siteofinterest,
  PFrecogsite,
  PRrecogsite,
  upstream = T,
  primersnumber = 5
)
}
\arguments{
\item{contig}{Sequence of contig that include site of interest}

\item{siteofinterest}{Approximate site for double strand break}

\item{PFrecogsite}{Sequence attached to forward strand of the target recognition seq.}

\item{PRrecogsite}{Sequence attached to reverse strand of the target recognition seq.}

\item{upstream}{Find site upstream or downstream}

\item{primersnumber}{How many primer pairs you want}
}
\value{
a list of primers, locations and sticky ends
}
\description{
The function finds the PAM sites and design primers to incorporate cleavage sites.
}
\details{
TBD
}
\examples{
CRISPRUSERprimercall(contig,site-50,PFrecogsite,PRrecogsite)

}
