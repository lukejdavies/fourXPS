EXP.4.Z<-function(z=z, mag=mag, band=band, col=col, mass=mass, sfr=sfr, agn=agn, res=0.25, passCut=-10.3, TSub=1200, MaxTime=6, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, offZ=0.05, reqProb=0.95, plot=T, ETCPath='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4most2/code/bin/', systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/', verbose=F) {
  
  id<-1
  SimName<-'EXP4Ztmp'
  make.Sim.Spec(id, z, mag, band, col, mass, sfr, agn, res=0.25, passCut=-10.3, simName=SimName)
  EXPSec<-TSub*c(1:ceiling(MaxTime*60*60/TSub))
  NSub<-c(1:ceiling(MaxTime*60*60/TSub))
  
  MulitSpec<-obs.Multi.4MOST(SimName, id, EXPSec, NSub, writeFITSOut=FALSE, airMass=airMass, IQ=IQ, skyBright=skyBright,tilt=tilt, misAlign=misAlign, ETCPath = ETCPath, systemModelDir=systemModelDir, verbose=verbose)
  
  system('rm -rf EXP4Ztmp')
  
  if (verbose==T){cat('Finished running ETC....' , '\n')}
  
  prob<-0
  off<-10
  for (i in 1:MulitSpec$Nconditions){
    if (off>offZ | prob<reqProb) {
  
      if (verbose==T){cat('Running auto.z with ',i,' of ',MulitSpec$Nconditions,' sub exposures' , '\n')}
        spec<-list(wave=MulitSpec$wave,flux=MulitSpec$conditions[[paste("condition",i,sep='')]]$flux, error=MulitSpec$conditions[[paste("condition",i,sep='')]]$error, xunit='ang', yunit='ang', z=MulitSpec$z)
        autozOut<-Autoz.Single.Spec(spec, verbose=F)
        spec$zFit<-autozOut$best$z
        spec$zProb<-autozOut$best$prob
        spec$zTemp<-autozOut$best$temp
        off<-abs(spec$zFit-MulitSpec$z)/MulitSpec$z
        prob<-spec$zProb
        NSub<-i
      }
    
  }

  if (off>offZ | prob<reqProb) {
    cat('****** WARNING NO GOOD SOLUTION FOUND ******', '\n')
    cat('Input Magnitude =', mag ,'\n')
    cat('Input Redshift =',spec$z ,'\n')
    cat('Fitted Redshift =',round(spec$zFit*1000)/1000 ,'\n')
    cat('Probability =',round(spec$zProb*1000)/1000 ,'\n')
    cat('# SubExposures =',NSub ,'\n')
    cat('Time for each SubExposure (min) =',TSub/60 ,'\n')
    cat('Total time (min) =',NSub*TSub/60  ,'\n')
    cat('Max time set to (min) =', MaxTime*60  ,'\n')
    cat('Sky Brightness =',skyBright  ,'\n')
    cat('Air Mass =', airMass  ,'\n')
    
  } else {
  
    if (verbose==T){
      cat('-----------------------' ,'\n')
      cat('GOOD REDSHIFT FOUND:' ,'\n')
      cat('-----------------------' ,'\n')
      cat('Input Magnitude =', mag ,'\n')
      cat('Input Redshift =',spec$z ,'\n')
      cat('Fitted Redshift =',round(spec$zFit*1000)/1000 ,'\n')
      cat('Probability =',round(spec$zProb*1000)/1000 ,'\n')
      cat('# SubExposures =',NSub ,'\n')
      cat('Time for each SubExposure (min) =',TSub/60 ,'\n')
      cat('Total time (min) =',round(NSub*TSub/60)  ,'\n')
      cat('Sky Brightness =',skyBright  ,'\n')
      cat('Air Mass =', airMass  ,'\n')
    }
  
  }
  
  if (plot==T){
    
    data('filtered_templates', package='auto.z')
    TempWave<-tdata$col[[5]][which(tdata$col[[1]] == autozOut$best$temp),]
    TempFlux<-tdata$col[[7]][which(tdata$col[[1]] == autozOut$best$temp),]
    
    TempFlux<-(TempFlux/max(TempFlux,na.rm=T))*(max(spec$flux,na.rm=T)/2.0)+median(spec$flux,na.rm=T)
    TempWave<-10^(TempWave)*(1+spec$zFit)
    
    spec.plot(spec, col='navy')
    lines(TempWave,TempFlux, col='indianred2')
    legend('topright', legend=c(paste('zIn = ',spec$z, sep=''),paste('zOut = ',round(spec$zFit*1000)/1000, sep=''),paste('Prob = ',round(spec$zProb*1000)/1000, sep=''),paste('Texp(min) = ',round(NSub*(TSub/60)),sep='')),bg='white')
    legend('topleft', legend=c('4MOST Spectrum', 'Template Fit'), text.col=c('navy','indianred2'),bg='white')
    
    if (off>offZ | prob<reqProb) {
      text(median(spec$wave,na.rm=T), median(spec$flux,na.rm=T)-(max(spec$flux,na.rm=T)-min(spec$flux,na.rm=T))/3.0, '** WARNING, NO GOOD SOLUTION FOUND **', cex=2, col='red')
    }
    
  }
  
  spec$TSub<-TSub
  spec$NSub<-NSub
  spec$EXPSec<-NSub*(TSub)
  
  
  spec$mag<-mag
  spec$band<-band
  spec$col<-col
  spec$mass<-mass
  spec$sfr<-sfr
  spec$agn<-agn
  spec$airMass<-airMass
  spec$IQ<-IQ
  spec$skyBright<-skyBright
  spec$tilt<-tilt
  spec$misAlign<-misAlign
  
  return=spec
    
}