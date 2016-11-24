\name{EXP.4.Z}
\alias{EXP.4.Z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find minimum exposure time required to obtain a redshift using simulated 4MOST observations.
}
\description{
High level function that combines, simulated spectrum generation (make.Sim.spec), mock 4MOST observations (using the 4MOST ETC) and redshifting (using auto.z). Code allows user to input simulated spectral paramters and 4MOST observing conditions, and will then return the minimum exposure time required to obtain a secure and correct redshift for the input paramters.  
}
\usage{
EXP.4.Z(z=z, mag=mag, band=band, col=col, mass=mass, sfr=sfr, agn=agn, res=0.25, passCut=-10.3, TSub=1200, MaxTime=6, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, offZ=0.05, reqProb=0.95, plot=T, ETCPath='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4most2/code/bin/', systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/', verbose=F) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
Redshift required for simulated spectrum 
}
  \item{mag}{
Fibre magnitude required for simulated spectrum 
}
  \item{band}{
Band in which magnitude is measured
}
  \item{col}{
g-i colour required for simulated spectrum 
}
  \item{mass}{
Log[M*/Msun] required for simulated spectrum 
}
  \item{sfr}{
SFR required for simulated spectrum 
}
  \item{agn}{
Identifier to assign AGN-like template to this source. F=do not use AGN template, N=use narrowline AGN template, B=use broadline AGN template.
}
  \item{res}{
Output resolution of the generated spectrum (Ang/pixel). Default=0.25
}
  \item{passCut}{
The specific star formation rate cut with which to use passive templates in Log[Msun/yr/Msun]. Default=-10.3. Thus input sources with log[sSFR] < -10.3 will use passive templates only.
}
  \item{TSub}{
Length of each SubExposure to observe, in seconds
}
  \item{MaxTime}{
Maximum time to observe, hours
}
  \item{airMass}{
Simulted Air mass
}
  \item{IQ}{
Simulted Image quality
}
  \item{skyBright}{
Simulted Sky Brightness (mags/arcsec^2)
}
  \item{tilt}{
Simulated Fibre Tilt (mm)
}
  \item{misAlign}{
Simulated fibre mis-allignment (mm)
}
  \item{offZ}{
Allowed offset between auto.z and input redshift for observation to be deemed sucessful
}
  \item{reqProb}{
Probability of auto.z redshift required for observation to be deemed sucessful
}


  \item{plot}{
Should resultant spectrum and best fit template be plotted? TRUE/FALSE
}
  \item{ETCPath}{
Location for 4FS-ETC executable
}
  \item{systemModelDir}{
Path to current 4MOST system model
}
  \item{verbose}{
Keep me informed of whats going on! TRUE/FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{flux }{4MOST observed spectral flux}
\item{wave }{4MOST observed Wavelength}
\item{error }{4MOST observed spectral error}
\item{xunit }{Units in the x-direction. Units are 'ang', 'hz', 'm', 'micron','nm'. }
\item{yunit }{Units in the y-direction. Units are 'ang'=ergs/sec/cm^2/Ang, 'hz'=ergs/sec/cm^2/Hz, 'Jy'=Jansky.}
\item{z }{Input redshift of simulated spectrum}
\item{zFit }{Best fit auto.z redshift}
\item{zProb }{Probability of best-fit redshift}
\item{zTemp }{Best fit auto.z template number}
\item{TSub }{Inout Subexposure length, sec}
\item{NSub }{Number of exposures to obtain redshift}
\item{EXPSec }{Total time of simualted observation}
\item{mag }{Input simulated fibre magnitude}
\item{band }{Band of simulated magnitude}
\item{col }{Simulated g-i colour}
\item{mass }{Simulated stellar mass  log[M*/M]}
\item{sfr }{Simulated SFRs}
\item{agn }{Simulated AGN identifier. F=non-AGN, N=use narrowline AGN template, B=use broadline AGN template.}
\item{airMass }{Air Mass of simulated observation}
\item{IQ }{Image quality of simulated observation}
\item{skyBright }{Sky Brightness of simulated observation}
\item{tilt }{Fibre tilt of simulated observation}
\item{misAlign }{Fibre Mis-alignment of simulated observation}

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

#This code will generate a single spectrum with redshift=0.3, fibre magnitude=22.5 in the 'r' band, g-i colour of 1.1, stellar mass=10^10 (M*/Msun) and SFR= 0.01 (i.e. this will use a passive template). The simulated spectrum will have an output resultion of 0.25 Ang/pix and does not use an AGN=like template. Observations are simulated with 1200sec subsexposires for a maximum obsering time of 6h. The maximum allowable offset between the auto.z-determined redshfit and input redshift is 0.01 and the auto-z redshift must be secure to 0.99 probability.   

Out<-EXP.4.Z(z=0.3, mag=22.0, band='r', col=1.1, mass=10.0, sfr=0.01, agn='F', res=0.25, TSub=1200, MaxTime=6, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, offZ=0.01, reqProb=0.99, plot=T, verbose=T)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }
\keyword{  }% __ONLY ONE__ keyword per line
