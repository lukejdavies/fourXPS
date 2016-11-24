Sim.to.4MOST<-function(z=z, mag=mag, band=band, col=col, mass=mass, sfr=sfr, agn=agn, res=0.25, passCut=-10.3, EXPSec=EXPSec, NSub=NSub, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, ETCPath='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4most2/code/bin/', systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/', verbose=F) { 
  
  id<-1
  outname='NA'
  simName='tmpSimTo4MOST'
  if (verbose==T){cat('Generating Spectrum..... ','\n')}
  
  make.Sim.Spec(id, z, mag, band, col, mass, sfr, agn, res=0.25, passCut=-10.3, simName=simName)
  
  spec<-obs.Single.4MOST(simName,id, EXPSec, NSub, OutName, writeFITSOut=FALSE, airMass=airMass, IQ=IQ, skyBright=skyBright,tilt=tilt, misAlign=misAlign, ETCPath=ETCPath,systemModelDir=systemModelDir, verbose=verbose)
  system('rm -rf tmpSimTo4MOST')
  
  spec[['col']]<-col
  spec[['mass']]<-mass
  spec[['sfr']]<-sfr
  spec[['agn']]<-agn
  
  return(spec)  
  
}