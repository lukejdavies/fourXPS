\name{make.ETCpar}
\alias{make.ETCpar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to generate 4MOST ETC input parmaters files
}
\description{
This function allows the easy generation of the two input parameter files needed for running the 4MOST 4FS-ETC. Some parameters are fixed as standard, but this code allows easy editing to run 4FS-ETC over a given simulated spectra (using simName+id), for a given set of observing conditions. The current 4MOST system model location must be provided using the 'systemModelDir' value.  
}
\usage{
make.ETCpar(simName, id, EXPSec, NSub, MagIn, MagSim, airMass = 1.4, IQ = 1.1, skyBright = 21.77, tilt = 6, misAlign = 0.1, systemModelDir = "/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simName}{
Simulation name of spectrum to define input files. Based on make.Sim.Spec outputs.
}
  \item{id}{
Spectrum ID within simulation to define input files. Based on make.Sim.Spec outputs.
}
  \item{EXPSec}{
Exposure time in seconds. Can be a scalar or vector
}
  \item{NSub}{
Number of sub exposures. Can be a scalar or vector, but must match EXPSec 
}
  \item{MagIn}{
Input fibre magnitude of simulated spectrum
}
  \item{MagSim}{
Fibre mangitudes to simulate. Can be a scalar value or a vector with the format: c(M1, M2, interval). This will produce input files to simulate the input source at a number of magnitude intervals between M1 and M2. 
}
  \item{airMass}{
Air mass to simulate. Can be a scalar or vector
}
  \item{IQ}{
Image quality to simulate (V-band,FWHM,arcsec). Can be a scalar or vector
}
  \item{skyBright}{
Sky brightness to simulate (mags/arcsec^2). Can be a scalar or vector
}
  \item{tilt}{
Fibre tilt to simulate (mm). Can be a scalar or vector
}
  \item{misAlign}{
Fibre misallignment to simulate (mm). Can be a scalar or vector
}

  \item{systemModelDir}{
Current 4MOST system model location
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

# Will generate a parameters file to simulate source ID=1 and simulation="test", with an input magnitude of 23, for 100sec with 1 subexposure. 
make.ETCpar('test',1, 100, 1, 23, 23, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/')

# Will generate a parameters file to simulate source ID=1 and simulation="test", with an input magnitude of 23. Simulations would be produced with the spectrum scaled to 21.1, 21.2.....21.9,22.0, each for 100sec with 1 subexposure, and 500sec with 5 subexposures. 
make.ETCpar('test',1, c(100,500), c(1,5), 23, c(21,22,0.1), airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/')

# Will generate a parameters file to simulate source ID=1 and simulation="test", with an input magnitude of 23, for 100sec with 1 subexposure. Simulations will be prodiced at three different airmasses and sky brightnesses. 
make.ETCpar('test',1, 100, 1, 23, 23, airMass=c(1.3,1.4,2.0), IQ=1.1, skyBright=c(20.5,21.77 22.2),tilt=6.0, misAlign=0.1, systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }
\keyword{  }% __ONLY ONE__ keyword per line
