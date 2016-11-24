\name{make.Sim.Spec}
\alias{make.Sim.Spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to generate mock simulated spectra
}
\description{
This function allows the user to generate multiple mock spectra at FITS files based on a series of redshift, magnitude, g-i colour, stellar mass, sfr and AGN-class inputs. 
}
\usage{
make.Sim.Spec(id, z, mag, band, col, mass, sfr, agn, res=0.25, passCut=-10.3, simName='testSim')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
Unique id for output spectra. Must be a vector with the same size as other main inputs (z, mag, band, col, mass, sfr, agn).
}
\item{z}{
 Unique redshift to simulate for each id.
 }
\item{mag}{
 Unique magnitude to simulate for each id.
 } 
\item{band}{
 Filter band associated with mag (u,g,r,i,z).
 }
 \item{col}{
 Unique rest-frame g-i colour to match to library of template spectra.
 }
 \item{mass}{
 Unique log stellar mass to match to library of template spectra (log[M*/Msun]).
 }
 \item{sfr}{
 Unique sfr to match to library of template spectra (Msun/yr).
 }
 \item{agn}{
Unique identifier to assign AGN-like template to this source. F=do not use AGN template, N=use narrowline AGN template, B=use broadline AGN template.
 }
 \item{res}{
Output resolution of the generated spectrum (Ang/pixel). Default=0.25
 }
 \item{passCut}{
The specific star formation rate cut with which to use passive templates in Log[Msun/yr/Msun]. Default=-10.3. Thus input sources with log[sSFR] < -10.3 will use passive templates only.  
 }
 \item{simName}{
 Output name of simulation. A directory will be created with this name and each spectrum will have this string in its filename. Default='testSim'
 }
 
 
 
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function outputs spectra with the form *outname*/*outName*_*id*.fits. One file for each unique id.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Luke Davies
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example generates 10 random spectra......

id<-seq(1:10) 

z<-runif(length(id),0.01,1) # Redshifts are randomly generated between 0.01 and 1

sfr<-runif(length(id),0.0,20) # SFRs are randomly generated between 0 and 20

sfr[round(runif(3,1,length(id)))]<-0.0 # 3 id's are foced to have SFR=0 (passive sytems)

mass<-runif(length(id),9,11) # Masses are randomly generated between log[M*/M]=9-11

mass[which(sfr==0)]<-runif(length(which(sfr==0)),10,11.5) # The passive systems forced above are also forced to have higher masses, log[M*/M]=10-11.5

col<-runif(length(id),0.3,1.05) # g-i colours are randomly generated between 0.3 and 1.05

mass[which(sfr==0)]<-runif(length(which(sfr==0)),1.1,1.25) # The passive systems forced above are also forced to have redder colours, 1.1-1.25

mag<-runif(length(id),20,23) # Magnitudes are randomly generated between 20 and 23

band<-'r'# The magnitude band is set to 'r'

agn<-id 

agn[]<-'F' # Initally all sources are set to not use AGN templates

agn[9]<-'B' # One sources is set to use a broad-line AGN template

agn[10]<-'N'# One sources is set to use a narrow-line AGN template


SimName<-'TestSim' #Putput name is set to 'TestSim'

make.Sim.Spec(id, z, mag, band, col, mass, sfr, agn, res=0.25, passCut=-10.3, simName=SimName)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
