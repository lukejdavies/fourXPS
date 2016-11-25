\name{Sim.to.4MOST}
\alias{Sim.to.4MOST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a single spectrum and mock observe using 4MOST
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Sim.to.4MOST(z=z, mag=mag, band=band, col=col, mass=mass, sfr=sfr, agn=agn, res=0.6, passCut=-10.3, EXPSec=EXPSec, NSub=NSub, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, ETCPath='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4most2/code/bin/', systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/', verbose=F)
  
}
%- maybe also 'usage' for other objects documented here.
\arguments{

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
\item{col }{Input spectrum g-i colour }
\item{mass }{Input spectrum log stellar mass }
\item{sfr }{Input spectrum sfr }
\item{agn }{Input spectrum AGN flag }
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

spec<-Sim.to.4MOST(z=0.1, mag=21, band='r', col=0.75, mass=10.5, sfr=5, agn='F', res = 0.6, passCut = -10.3, EXPSec=9000, NSub=8, airMass = 1.4, IQ = 1.1, skyBright = 21.77, tilt = 6, misAlign = 0.1,verbose = F)
spec.plot(spec)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }
\keyword{  }% __ONLY ONE__ keyword per line
