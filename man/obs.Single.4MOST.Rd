\name{obs.Single.4MOST}
\alias{obs.Single.4MOST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for simulating a single spectrum observation with 4MOST using outputs from make.Sim.Spec. 
}
\description{
This function takes a single spectrum output from make.Sim.Spec (using simName and id), and passes the spectrum thought the 4MOST exposure time calculator for a set of desired observing conditions. The function requires an installed version of the 4MOST ETC, which runs from the command line with using '4FS_ETC'. For 4MOST consortium members, the 4MOST ETC and current system design is availble here: http://wiki.4most.eu/4most-facility-simulator#toc10 .  

}
\usage{
obs.Single.4MOST(simName,id, EXPSec, NSub, OutName, writeFITSOut=TRUE, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, ETCPath='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4most2/code/bin/', systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/', verbose=F)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simName}{
The simulation name defined in make.Sim.Spec.
}
  \item{id}{
The id of the spectrum to simulate, also defined in make.Sim.Spec. 
}
  \item{EXPSec}{
Simulated exposure time in seconds. 
}
  \item{NSub}{
Number of subexposures in EXPSec.
}
  \item{OutName}{
Output string value name of simulated spectrum
}
  \item{writeFITSOut}{
Should output spectrum be saved as FITS? TRUE/FALSE. This produces a FITs file with the simulated spectra. Within this file row1=4MOST simulated flux, row2=4MOST simulated error, and row3=input flux (at top of atmosphere).
}
  \item{airMass}{
Simulated airmass.
}
  \item{IQ}{
Simulated image quality (V-band,FWHM,arcsec)
}
  \item{skyBright}{
Simulated sky brightness (mag/arcsec^2)
}
  \item{tilt}{
Simulated fibre tilt (mm)
}
  \item{misAlign}{
Simulated fibre misallignment (mm)
}

  \item{ETCPath}{
Path the ETC executable
}

  \item{systemModelDir}{
String to directory path of current 4MOST system model (provided as part of 4FS)
}
  \item{verbose}{
Should the ETC provide verbose comments? TRUE/FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{wave }{Wavelength of simulated spectrum}
\item{flux }{Flux of simulated spectrum}
\item{error }{Error of simulated spectrum}
\item{inputFlux }{Input flux of simulated spectrum (top of atmosphere)}
\item{id }{ID of simulated spectrum}
\item{sim }{Simulation name}
\item{z }{Redshift}
\item{band }{Band in which magnitude is defined}
\item{mag }{Mangitude}
\item{xunit }{Units in the x-direction. Allowed units are 'ang', 'hz', 'm', 'micron','nm'. }
\item{yunit }{Units in the y-direction. Allowed units are 'ang'=ergs/sec/cm^2/Ang, 'hz'=ergs/sec/cm^2/Hz, 'Jy'=Jansky.}
\item{EXPSec }{Exposure time in seconds}
\item{NSub }{Number of sub exposures}
\item{axDat }{The axDat information used in writeFITS - see this package for details}
\item{airmass }{Input 4MOST airmass }
\item{IQ }{Input 4MOST image quality}
\item{skyBright }{Input 4MOST skybrightness }
\item{tilt }{Input 4MOST fibre tilt }
\item{misAlign }{Input 4MOST fibre misallignment }

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

spec<-obs.Single.4MOST("TestSim", 1, 4800,4,'test.fits', writeFITSOut=FALSE, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/', verbose=F)

# Will simulate spectrum ID=1 for simulation, "TestSim". Will expose for 4800 seconds in 4 subexposures (each of 1200sec). Will also save the spectrum as "test.fits".

spec.plot(spec)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
\keyword{  }% __ONLY ONE__ keyword per line
