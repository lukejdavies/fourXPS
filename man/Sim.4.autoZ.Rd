\name{Sim.4.autoZ}
\alias{Sim.4.autoZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produce a single simulated spectrum, mock observe with 4MOST, and determine automatic redshift fit.
}
\description{
Function for automatically mock observing (using 4MOST) and redshifting a single simulated spectrum. Allows the uses to test a single spectrum type and magnitude to see if it will be sucessfully observed by 4MOST in a given set of observing conditions. 
}
\usage{
Sim.4.autoZ(z = z, mag = mag, band = band, col = col, mass = mass, sfr = sfr, agn = agn, res = 0.25, passCut = -10.3, EXPSec = EXPSec, NSub = NSub, airMass = 1.4, IQ = 1.1, skyBright = 21.77, tilt = 6, misAlign = 0.1, plot = T, ETCPath = "/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4most2/code/bin/", systemModelDir = "/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/", verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
Redshift required for simulated spectrum 
}
  \item{mag}{
Magnitude required for simulated spectrum 
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
  \item{EXPSec}{
Total exposure time to Simulate in sec
}
  \item{NSub}{
Number of sub exposures in EXPSec 
}
  \item{airMass}{
Simulted Air mass
}
  \item{IQ}{
Simulted Image quality
}
  \item{skyBright}{
Simulted Sky Brightness
}
  \item{tilt}{
Simulated Fibre Tilt
}
  \item{misAlign}{
Simulated fibre mis-allignment
}
  \item{plot}{
Should resultant spectrum and fit be plotted? TRUE/FALSE
}
  \item{ETCPath}{
Location for $FS-ETC executable
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
\item{zFit }{Best Fit auto.z  redshift}
\item{zProb }{Probability of best fit auto.z redshift}
\item{zTemp }{template number of best fit auto.z redshift}

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

spec<-Sim.4.autoZ(z=0.1, mag=21, band='r', col=0.75, mass=10.5, sfr=5, agn='F', res = 0.25, passCut = -10.3, EXPSec=9000, NSub=8, airMass = 1.4, IQ = 1.1, skyBright = 21.77, tilt = 6, misAlign = 0.1, plot=T, verbose = F)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }
\keyword{  }% __ONLY ONE__ keyword per line
