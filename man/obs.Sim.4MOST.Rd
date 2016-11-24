\name{obs.Sim.4MOST}
\alias{obs.Sim.4MOST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mock observe a set of simulated sepectra using the 4MOST ETC.  
}
\description{
Function allows a simulated sample (produced by make.Sim.Spec) to be mock observed by 4MOST. Either a single set of observing spec can be applied to all sources, or vectors provided for each individual source.
}
\usage{
obs.Sim.4MOST(simName, id, EXPSec, NSub, airMass = 1.4, IQ = 1.1, skyBright = 21.77, tilt = 6, misAlign = 0.1,  writeFITSOut=FALSE, ETCPath = "/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4most2/code/bin/", systemModelDir = "/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/", verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simName}{
Input name of simulation as produced in make.Sim.Spec
}
  \item{id}{
Scalar or Vector of IDs of sepctrum in simulation as produced in make.Sim.Spec
}
  \item{EXPSec}{
Exposure times to simulate in sec. Either scalar to use the same exposure time for all sources, or a vector of the same length as id.
}
  \item{NSub}{
Number of subexposures within EXPSec. Either scalar to use the same NSub for all sources, or a vector of the same length as id.
}
  \item{airMass}{
Air Masses to simulate. Scalar or verctor with length of id. 
}
  \item{IQ}{
Image quality to simulate. Scalar or verctor with length of id. 
}
  \item{skyBright}{
Sky Brightnesses to simulate (mags/arcsec^2). Scalar or verctor with length of id. 
}
  \item{tilt}{
Fibre tilts to simulate (mm). Scalar or verctor with length of id. 
}
  \item{misAlign}{
Fibre mis-allignaments to simulate (mm). Scalar or verctor with length of id. 
}
  \item{writeFITSOut}{
Should resultant spectra be written as FITS. TRUE/FALSE
}
  \item{ETCPath}{
Location of 4FS-ETC excecuatable
}
  \item{systemModelDir}{
Path to curretn 4MOST system model
}
  \item{verbose}{
Keep me informed of whats going on! TRUE/FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{Sim }{Input simulation name}
  \item{NSpec }{Number of spectra simulated}
  \item{specX }{Nested spectrum-specific outputs (X is spectrum number...)}
  \item{specX$id }{Input id of spectrum}
   \item{specX$flux }{Output 4MOST observed flux for given spectrum}
   \item{specX$error }{Output 4MOST observed error for given spectrum}
   \item{specX$ExpSec }{Output 4MOST exposure time for given spectrum, in seconds}
   \item{specX$NSub }{Output 4MOST Number of sub exposures for given spectrum}
   \item{specX$airmass }{Output 4MOST airmass for given spectrum}
   \item{specX$IQ }{Output 4MOST image quality for given spectrum}
   \item{specX$skyBright }{Output 4MOST skybrightness for given spectrum}
   \item{specX$tilt }{Output 4MOST fibre tilt for given spectrum}
   \item{specX$misAlign }{Output 4MOST fibre misallignment for given spectrum}
   \item{specX$FluxOrig2 }{Input flux from spectrum, interpolated onto output wavelength scale.}
   \item{specX$z }{Redshift of input spectrum}
   \item{specX$band }{Band of input magnitude}
   \item{specX$mag }{Input magnitude}
   \item{specX$xunit }{xunit of spectrum, see spec.tools/get.spec}
   \item{specX$yunit }{yunit of spectrum, see spec.tools/get.spec}
   \item{specX$axDat }{axDat information for current spectrum, see FITSio/readFITS}



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

simName<-'TestSim'
id<-c(1:10)
NSub<-round(runif(10,1,10))
EXPSec<-1200*NSub
out<-obs.Sim.4MOST(simName=simName,id=id, EXPSec=EXPSec, NSub=NSub, airMass = 1.4, IQ = 1.1, skyBright = 21.77, tilt = 6, misAlign = 0.1, writeFITSOut=FALSE, ETCPath = "/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4most2/code/bin/", systemModelDir = "/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/", verbose = F)

cols<-rainbow(out$NSpec)
for (i in 1:out$NSpec){
  spec<-list(wave=out[[paste("spec",i,sep='')]]$wave,flux=out[[paste("spec",i,sep='')]]$flux, xunit='ang', yunit='ang', z=out$z)
  if (i == 1 ){spec.plot(spec, lwd=0.5,col=cols[i])}
  if (i > 1 ){spec.plot(spec, lwd=0.5,col=cols[i],oplot=T)}
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }
\keyword{  }% __ONLY ONE__ keyword per line
