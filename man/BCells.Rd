\name{BCells}
\alias{BCells}
\docType{data}
\title{
Two-Photon Data: B Cells in a Lymph Node
}
\description{
Labelled B cells were adoptively transfered and intravitally
imaged (using two-photon microscopy) inside a peripheral lymph
node of the recipient mouse. These data illustrate the characteristic
"random-walk-like" motion pattern of B cells in lymph nodes.
}
\usage{data("BCells")}
\format{
  An S3 object of class "tracks"; a list with 24 elements. Each 
  element name identifies a cell track. Each element is a matrix
  containing the following four columns.
  \describe{
    \item{\code{t}}{the time (in seconds)}
    \item{\code{x}}{The X coordinate (in micrometers)}
    \item{\code{y}}{The Y coordinate (in micrometers)}
    \item{\code{z}}{The Z coordinate (in micrometers)}
  }
}
\source{
Data were generated in 2012 in the Mark J. Miller Lab, 
Department of Medicine,
Washington University in St Louis, USA.
}
\references{
Zinselmeyer BH, Dempster J, Wokosin DL, Cannon JJ, Pless R, Parker I and Miller MJ
(2009), 
Two-photon microscopy and multi-dimensional analysis of cell dynamics.
\emph{Methods in Enzymology}, \bold{461}:349--78.
doi:10.1016/S0076-6879(09)05416-0

Konjufca V and Miller MJ (2009),
Imaging Listeria monocytogenes infection in vivo. 
\emph{Current Topics in Microbiology and Immunology}, \bold{334}:199--226.   
doi:10.1007/978-3-540-93864-4_9

Kreisel D, Nava RG, Li W, Zinselmeyer BH, Wang B, Lai J, Pless R, Gelman AE, Krupnick AS, 
and Miller MJ (2010),
In vivo two-photon imaging reveals monocyte-dependent neutrophil extravasation during 
pulmonary inflammation. 
\emph{PNAS}, \bold{107}(42):18073--18078.
doi:10.1073/pnas.1008737107
}
\examples{
## load the tracks
data(BCells)
## visualize the tracks (calls function plot.tracks)
plot(BCells)
}
\keyword{datasets}
