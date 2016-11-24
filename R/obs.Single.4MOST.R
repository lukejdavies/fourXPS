obs.Single.4MOST<-function(simName,id, EXPSec, NSub, OutName, writeFITSOut=TRUE, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, ETCPath='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4most2/code/bin/', systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/', verbose=F){

      
      system(paste('mkdir tmpETC',simName, sep=''))  
  
      red_thru<-readFITS(paste(systemModelDir,'LRS/lrs_red_material_4fs_efficiency_total.fits', sep=''), hdu=1)
      blue_thru<-readFITS(paste(systemModelDir,'LRS/lrs_blue_material_4fs_efficiency_total.fits', sep=''), hdu=1)
      green_thru<-readFITS(paste(systemModelDir,'LRS/lrs_green_material_4fs_efficiency_total.fits', sep=''), hdu=1)
  
   
          spec<-readFITS(paste(simName,'/',simName,'_',id,'.fits', sep=''))
          FluxOrig<-spec$imDat
          WaveOrig<-spec$axDat$crval[1]+((seq(1,length(spec$imDat),1)-spec$axDat$crpix[1])*spec$axDat$cdelt[1])
          specMag<-as.numeric(spec$hdr[which(spec$hdr=='MAG')+1])
          specZ<-as.numeric(spec$hdr[which(spec$hdr=='Z')+1])
          band<-spec$hdr[which(spec$hdr=='BAND')+1]
          specT<-list(flux=FluxOrig, wave=WaveOrig, xunit='ang', yunit='ang')
          MagIN<-magABspec(specT, filter='r')
          
          
          axDat<-spec$axDat
          
          system(paste('rm -rf tmpETC',simName,'/specout*.fits',sep=''))
          
          if (verbose==T){cat('Producing ETC input files..... ','\n')}
          
          make.ETCpar(simName,id, EXPSec, NSub, MagIN, MagIN, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, systemModelDir=systemModelDir)
            
          if (verbose==T){cat('Running ETC..... ','\n')}
          
          if (verbose==F) {system(paste(ETCPath,'4FS_ETC PARAM_FILENAME=',simName,'_ETC_input_params_tmp.txt \ >& ',simName,'_etc.log &',sep=''))}
          if (verbose==T) {system(paste(ETCPath,'4FS_ETC PARAM_FILENAME=',simName,'_ETC_input_params_tmp.txt &',sep=''))}
          Sys.sleep(2)
  
          conditions<-readFITS(paste('tmpETC',simName,'/specout_template_',simName,'_',id,'_LRS_red.fits',sep=''), hdu=1)
          results_r<-readFITS(paste('tmpETC',simName,'/specout_template_',simName,'_',id,'_LRS_red.fits',sep=''), hdu=2)
          results_b<-readFITS(paste('tmpETC',simName,'/specout_template_',simName,'_',id,'_LRS_blue.fits',sep=''), hdu=2)
          results_g<-readFITS(paste('tmpETC',simName,'/specout_template_',simName,'_',id,'_LRS_green.fits',sep=''), hdu=2)

          
          system(paste('rm -rf tmpETC',simName, sep=''))
          system(paste('rm -rf ',simName,'_ETC_input_params_tmp.txt',sep=''))
          system(paste('rm -rf ',simName,'_tmp.ETC', sep=''))
          system(paste('rm -rf ',simName,'_etc.log', sep=''))
          
          if (verbose==T){cat('Combining spectrum arms...... ','\n')}
          
          red_thru2<-approx(red_thru$col[[1]]*10, red_thru$col[[12]], results_r$col[[1]])$y
          blue_thru2<-approx(blue_thru$col[[1]]*10, blue_thru$col[[12]], results_b$col[[1]])$y
          green_thru2<-approx(green_thru$col[[1]]*10, green_thru$col[[12]], results_g$col[[1]])$y
  
          red_flux<-(results_r$col[[3]])
          green_flux<-(results_g$col[[3]])
          blue_flux<-(results_b$col[[3]])
  
          
          red_wave<-results_r$col[[1]]
          green_wave<-results_g$col[[1]]
          blue_wave<-results_b$col[[1]]
          
          
          red_noise<-results_r$col[[4]]
          green_noise<-results_g$col[[4]]
          blue_noise<-results_b$col[[4]]
  
          red_SNR<-results_r$col[[5]]
          green_SNR<-results_g$col[[5]]
          blue_SNR<-results_b$col[[5]]
          
          
          
          red_error<-red_flux
          green_error<-green_flux
          blue_error<-blue_flux
  
          red_flux2<-red_flux
          green_flux2<-green_flux
          blue_flux2<-blue_flux
          
          
          for (j in 1:length(red_flux)){
              red_error[j]<-rnorm(1,mean=0, sd=1.0*red_noise[j])
              red_flux2[j]<-red_flux[j]+red_error[j]                                      
          }
          for (j in 1:length(green_flux)){
              green_error[j]<-rnorm(1,mean=0, sd=1.0*green_noise[j])
              green_flux2[j]<-green_flux[j]+green_error[j]
  
          }
          for (j in 1:length(blue_flux)){
              blue_error[j]<-rnorm(1,mean=0, sd=1.0*blue_noise[j])
              blue_flux2[j]<-blue_flux[j]+blue_error[j]
          }
  
          red_flux2<-red_flux2/red_thru2
          green_flux2<- green_flux2/green_thru2
          blue_flux2<- blue_flux2/blue_thru2
          
  
          red_flux2<-red_flux2[350:length(red_wave)]
          red_error<-red_error[350:length(red_wave)]
          red_wave<-red_wave[350:length(red_wave)]
  
          green_flux2<-green_flux2[350:(length(green_wave)-350)]
          green_error<-green_error[350:(length(green_wave)-350)]
          green_wave<-green_wave[350:(length(green_wave)-350)]
  
          blue_flux2<-blue_flux2[0:(length(blue_wave)-350)]
          blue_error<-blue_error[0:(length(blue_wave)-350)]
          blue_wave<-blue_wave[0:(length(blue_wave)-350)]
  
  
          
          redFlux_sel1<-red_flux2[which(red_wave<=max(green_wave))]
          greenFlux_sel1<-green_flux2[which(green_wave>=min(red_wave))]
          tmp<-approx(red_wave[which(red_wave<=max(green_wave))],redFlux_sel1, green_wave[which(green_wave>=min(red_wave))])$y
          green_sc<-median(tmp-greenFlux_sel1, na.rm=T)
          green_flux2<-green_flux2+green_sc
          
          
          greenFlux_sel1<-green_flux2[which(green_wave<=max(blue_wave))]
          blueFlux_sel1<-blue_flux2[which(blue_wave>=min(green_wave))]
          tmp<-approx(green_wave[which(green_wave<=max(blue_wave))],greenFlux_sel1, blue_wave[which(blue_wave>=min(green_wave))])$y
          blue_sc<-median(tmp-blueFlux_sel1, na.rm=T)
          blue_flux2<-blue_flux2+blue_sc
  
          red_flux2<-red_flux2[100:length(red_wave)]
          red_error<-red_error[100:length(red_wave)]
          red_wave<-red_wave[100:length(red_wave)]
  
          green_flux2<-green_flux2[150:(length(green_wave)-100)]
          green_error<-green_error[150:(length(green_wave)-100)]
          green_wave<-green_wave[150:(length(green_wave)-100)]
  
          blue_flux2<-blue_flux2[0:(length(blue_wave)-150)]
          blue_error<-blue_error[0:(length(blue_wave)-150)]
          blue_wave<-blue_wave[0:(length(blue_wave)-150)]
          

          spec<-list()
          spec$wave<-c(blue_wave, green_wave,red_wave)
          spec$flux<-c(blue_flux2, green_flux2,red_flux2)
          spec$error<-abs(c(blue_error, green_error,red_error))  
          
          res<-median(c(spec$wave,0)-c(0,spec$wave),na.rm=T)
          wave<-seq(min(spec$wave),max(spec$wave),res)
          flux<-approx(spec$wave,spec$flux,wave)$y
          error<-approx(spec$wave,spec$error,wave)$y
          FluxOrig<-approx(WaveOrig,FluxOrig,wave)$y
          
  
          if (writeFITSOut==TRUE){
            
            if (verbose==T){cat('Writing FITS files...... ','\n')}
            
            axDat$crval[1]<-min(wave)
            axDat$crpix[1]<-1
            axDat$cdelt[1]<-res
            axDat$ctype[1]<-'LINEAR'
            axDat$cunit[1]<-'Angstrom'
            axDat$len[1]<-length(wave)
            axDat$len[2]<-3
      
  
            writeFITSim(cbind(flux,error,FluxOrig), file = OutName, axDat = axDat, type='single')
            write.fitskey('OBJECT',id, OutName, comment = 'Input id', hdu = 1)    
            write.fitskey('SIM',simName, OutName, comment = 'Input Simulation identifier', hdu = 1)  
            write.fitskey('Z',specZ, OutName, comment = 'redshift', hdu = 1)
            write.fitskey('BAND',band, OutName, comment = 'band', hdu = 1)
            write.fitskey('MAG',specMag, OutName, comment = 'magnitude in band', hdu = 1)
            write.fitskey('CRTYPE1','WAVE', OutName, comment = 'X-axis type', hdu = 1)
            write.fitskey('AUNIT','Angstrom', OutName, comment = 'wavelength unit', hdu = 1)
            write.fitskey('BUNIT','erg/s/cm2/Angstrom', OutName, comment = 'flux density unit', hdu = 1)
            write.fitskey('EXP',EXPSec, OutName, comment = 'Exposure Time, Sec', hdu = 1)
            write.fitskey('NSUB',NSub, OutName, comment = 'Number of Sub Exposures', hdu = 1)
            write.fitskey('AIRMASS',airMass, OutName, comment = 'Simulated AirMass', hdu = 1)
            write.fitskey('IQ',IQ, OutName, comment = 'Simulated Image Quality', hdu = 1)
            write.fitskey('SKYBRI',skyBright, OutName, comment = 'Simulated Sky Brightness', hdu = 1)
            write.fitskey('TILT',tilt, OutName, comment = 'Simulated Fibre Tilt', hdu = 1)
            write.fitskey('MISALIG',misAlign, OutName, comment = 'Simulated Fibre Misalignment', hdu = 1)
            
           
          }
          

          spec<-list(wave=wave,flux=flux,error=error, inputFlux=FluxOrig, id=id, sim=simName, z=specZ,band=band,mag=specMag, xunit='ang',yunit='ang', EXPSec=EXPSec, NSub=NSub, axDat=axDat, airMass=airMass, IQ=IQ, skyBright=skyBright,tilt=tilt, misAlign=misAlign)
          return(spec)
}
