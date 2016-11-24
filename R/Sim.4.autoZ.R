Sim.4.autoZ<-function(z=z, mag=mag, band=band, col=col, mass=mass, sfr=sfr, agn=agn, res=0.25, passCut=-10.3, EXPSec=EXPSec, NSub=NSub, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, plot=T, ETCPath='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4most2/code/bin/', systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/', verbose=F) {
  
  
  spec<-Sim.to.4MOST(z=z, mag=mag, band=band, col=col, mass=mass, sfr=sfr, agn=agn, res = res, passCut = -passCut, EXPSec=EXPSec, NSub=NSub, airMass = airMass, IQ = IQ, skyBright = skyBright, tilt = tilt, misAlign = misAlign, ETCPath=ETCPath, systemModelDir = systemModelDir, verbose = verbose)
  
  if (verbose==T){cat('Running auto.z.....', '\n')}
  
  autozOut<-Autoz.Single.Spec(spec, verbose=F)

  spec$zFit<-autozOut$best$z
  spec$zProb<-autozOut$best$prob
  spec$zTemp<-autozOut$best$temp
  
  if (plot==T) {
    data('filtered_templates', package='auto.z')
    TempWave<-tdata$col[[5]][which(tdata$col[[1]] == autozOut$best$temp),]
    TempFlux<-tdata$col[[7]][which(tdata$col[[1]] == autozOut$best$temp),]
  
    TempFlux<-(TempFlux/max(TempFlux,na.rm=T))*(max(spec$flux,na.rm=T)/2.0)+median(spec$flux,na.rm=T)
    TempWave<-10^(TempWave)*(1+spec$zFit)
    
    spec.plot(spec, col='navy')
    lines(TempWave,TempFlux, col='indianred2')
    legend('topright', legend=c(paste('zIn = ',spec$z, sep=''),paste('zOut = ',round(spec$zFit*1000)/1000, sep=''),paste('Prob = ',round(spec$zProb*1000)/1000, sep=''),'4MOST Spectrum', 'Template Fit'), text.col=c('black', 'black', 'black','navy','indianred2'))
  }
  
  return=spec
  
}