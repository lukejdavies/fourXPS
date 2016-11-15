make.Sim.Spec<-function(id, z, mag, band, col, mass, sfr, agn, res=0.25, passCut=-10.3, outName='testSim'){

    system(paste('mkdir ', outName, sep=''))
    load(paste(.libPaths(),'/fourXPS/data/modelSpecTab.Rdata',sep=''))
    logsSFR<-log10(modelTab$SFR/(10.^modelTab$LOGMSTAR))
    HA_EW<-modelTab$HA_EW
    spec_a<-readFITS(paste(.libPaths(),'/fourXPS/data/',modelTab$MODNAME[1],sep=''))
    axDat<-spec_a$axDat

    sSFR<-sfr/(10.^mass)
    
    spm<-read.table(paste(.libPaths(),'/fourXPS/data/csp_Zeq1p0_chab.sed',sep=''), header=F)
    age<-unique(spm[,1])
    spm2<-read.table(paste(.libPaths(),'/fourXPS/data/csp_Zeq0p2_chab.sed',sep=''), header=F)
    age2<-unique(spm2[,1])

    tmpArr<-array(NA,dim=c(6252,2,(length(age)+length(age2))))

    count<-1
    
    for (i in 1:ceiling((length(tmpArr[1,1,])/2.0))){

        spm_wave<-spm[which(spm[,1]==age[i]),2]    
        spm_flux<-spm[which(spm[,1]==age[i]),3]
        spm_flux<-spm_flux[which(spm_wave>800 & spm_wave<10100)]
        spm_wave<-spm_wave[which(spm_wave>800 & spm_wave<10100)]
        
        tmpArr[,1,count]<-spm_wave
        tmpArr[,2,count]<-spm_flux

        count<-count+1

        spm_wave<-spm2[which(spm2[,1]==age2[i]),2]    
        spm_flux<-spm2[which(spm2[,1]==age2[i]),3]
        spm_flux<-spm_flux[which(spm_wave>800 & spm_wave<10100)]
        spm_wave<-spm_wave[which(spm_wave>800 & spm_wave<10100)]
        
        tmpArr[,1,count]<-spm_wave
        tmpArr[,2,count]<-spm_flux
        
        count<-count+1

    }



    
    for (i in 1:length(id)){

        

        if (agn[i]=='F'){

            #### Match Sim to Mod ####
            dist_col<-abs(col[i]-modelTab$COL)/sd(modelTab$COL, na.rm=T)
            dist_mass<-abs(mass[i]-modelTab$LOGMSTAR)/sd(modelTab$LOGMSTAR[which(is.finite(modelTab$LOGMSTAR)==T)], na.rm=T) 
            dist_sfr<-abs(log10(sfr[i])-log10(modelTab$SFR))/sd(log10(modelTab$SFR), na.rm=T)


            
            if (log10(sSFR[i]) > passCut) {
                dist_all<-sqrt(dist_col^2+dist_mass^2+dist_sfr^2)
                match<-which(dist_all==min(dist_all[which(logsSFR > -10.0 & HA_EW>=0.1)], na.rm=T))
            }
            if (log10(sSFR[i]) <= passCut) {
                dist_all<-sqrt(dist_col^2+dist_mass^2)
                match<-which(dist_all==min(dist_all[which(logsSFR < -10.5 & HA_EW<0.1)]))
            }
            ##########################
            

            ### Get correct model ####       
            spec<-get.spec(paste(.libPaths(),'/fourXPS/data/',modelTab$MODNAME[match],sep=''), xunit='ang', yunit='ang', z=0)
        }

        if (agn[i]!='F'){
            agn_reds<-c(0.1,0.3,0.5,0.7,1.1,1.3,1.55,1.85,2.25,2.75,3.5,4.5,6.0)
            agn_redsN<-c('010','030','050','070','110','130','155','185','225','275','350','450','600')
            z_match<-agn_redsN[which(abs(z[i]-agn_reds)==min(abs(z[i]-agn_reds)))]

            if (agn[i]=='B'){name<-paste(.libPaths(),'/fourXPS/data/AGN_v1.0_z',z_match,'_mr165_type1.fits',sep='')}
            if (agn[i]=='N'){name<-paste(.libPaths(),'/fourXPS/data/AGN_v1.0_z',z_match,'_mr165_type2.fits',sep='')}       
            spec<-get.spec(name, xunit='ang', yunit='ang', z=agn_reds[which(abs(z[i]-agn_reds)==min(abs(z[i]-agn_reds)))])
            spec$wave<-spec$wave/(1+spec$z)
            spec$z<-0
        }

        ### Scale spec in z ###
        specNew<-spec
        specNew$wave<-spec$wave*(1+z[i])
        specNew$z<-z[i]

        ### Interpolate back to sensible wavelength range ###
        specNew$flux<-approx(specNew$wave,specNew$flux,seq(3000,10000,res))$y
        specNew$wave<-seq(3000,10000,res)

        specNew$flux[which(is.finite(specNew$flux)==F)]<-0.0

        axDat$crval[1]<-3000
        axDat$cdelt[1]<-res
        axDat$ctype[1]<-'LINEAR'
        axDat$cunit[1]<-'Angstrom'
        axDat$len[1]<-length(specNew$flux)
        axDat$len[2]<-1
        
        ### Scale to desired r-mag ###
        wavefac=1e-10
        c<-299792458

        if (z[i]<0.6) {filter=getfilt('r')[,2:3]}
        if (z[i]>=0.6) {filter=getfilt('i')[,2:3]}
        filt_wave <- filter[,1]
        filt_trans <- filter[,2]

        interp<-approx(filt_wave, filt_trans, specNew$wave)
        filt_trans_interp <- interp$y
        filt_trans_interp[is.na(filt_trans_interp)] = 0

        fluxnu=(wavefac*specNew$flux*specNew$wave^2)/c

        flux_conv <- fluxnu*filt_trans_interp

        temp <- fluxnu*filt_trans_interp*specNew$wave
        temp2 <- filt_trans_interp*specNew$wave
        flux_filt <- sum(temp[which(is.finite(temp)==T)])/sum(temp2[which(is.finite(temp2)==T)])

        flux_need <- 10.^(-0.4*(mag[i]+48.6))

        flux_sc_Hz <- fluxnu*(flux_need/flux_filt)
        flux_sc <- (flux_sc_Hz*(2.998e18))/(specNew$wave^2)
       

        
       ##########################################
        

        chi_min<-1e9
        
        for (k in 200:length(tmpArr[1,1,])){
            spm_wave<-tmpArr[,1,k]
            spm_flux<-tmpArr[,2,k]
            spm_wave<-spm_wave*(1+z[i])

            if (z[i]<0.6) {filter=getfilt('r')[,2:3]}
            if (z[i]>=0.6) {filter=getfilt('i')[,2:3]}
            
           
            filt_wave <- filter[,1]
            filt_trans <- filter[,2]
            
            interp<-approx(filt_wave, filt_trans, spm_wave)
            filt_trans_interp <- interp$y
            filt_trans_interp[is.na(filt_trans_interp)] = 0

            fluxnu=(wavefac*spm_flux*spm_wave^2)/c

            flux_conv <- fluxnu*filt_trans_interp

            temp <- fluxnu*filt_trans_interp*spm_wave
            temp2 <- filt_trans_interp*spm_wave
            flux_filt <- sum(temp[which(is.finite(temp)==T)])/sum(temp2[which(is.finite(temp2)==T)])
            flux_need <- 10.^(-0.4*(mag[i]+48.6))

            flux_sc_Hz <- fluxnu*(flux_need/flux_filt)
            spm_flux_sc <- (flux_sc_Hz*(2.998e18))/(spm_wave^2)

        
            sel<-which(flux_sc>0)
            flux_sc_interp<-approx(spm_wave,spm_flux_sc,specNew$wave[sel])$y
            chi<-sum((flux_sc_interp-flux_sc[sel])^2/flux_sc[sel])

            if (chi<chi_min){
                sel2<-k
                chi_min<-chi
            }
            
         
        }

        specComp<-list(wave=specNew$wave,flux=flux_sc)
        i_mag<-magABspec(specComp,filter='i')
        
        spm_wave<-tmpArr[,1,sel2]
        spm_flux<-tmpArr[,2,sel2]
        spm_wave<-spm_wave*(1+z[i])

        filter=getfilt('i')[,2:3]
        filt_wave <- filter[,1]
        filt_trans <- filter[,2]

        interp<-approx(filt_wave, filt_trans, spm_wave)
        filt_trans_interp <- interp$y
        filt_trans_interp[is.na(filt_trans_interp)] = 0

        fluxnu=(wavefac*spm_flux*spm_wave^2)/c

        flux_conv <- fluxnu*filt_trans_interp

        temp <- fluxnu*filt_trans_interp*spm_wave
        temp2 <- filt_trans_interp*spm_wave
        flux_filt <- sum(temp[which(is.finite(temp)==T)])/sum(temp2[which(is.finite(temp2)==T)])
        flux_need <- 10.^(-0.4*(i_mag+48.6))

        flux_sc_Hz <- fluxnu*(flux_need/flux_filt)
        spm_flux_sc <- (flux_sc_Hz*(2.998e18))/(spm_wave^2)


        if (length(which(flux_sc==0 & specNew$wave<specNew$wave[which(flux_sc==max(flux_sc,na.rm=T))])>0)) {
            lowCut<-max(specNew$wave[which(flux_sc==0 & specNew$wave<specNew$wave[which(flux_sc==max(flux_sc,na.rm=T))])])
            low_sel<-which(specNew$wave<=(lowCut+500))
            spm_low_sel<-which(spm_wave<=(lowCut+500))
            spm_low_interp<-approx(spm_wave,spm_flux_sc,specNew$wave[low_sel])$y
            off_low<-median(spm_low_interp[which(flux_sc[low_sel]>0)]-flux_sc[which(flux_sc[low_sel]>0)],na.rm=T)
            spm_low_interp<-spm_low_interp-off_low
            low_sel2<-which(specNew$wave<=(lowCut+10))
            flux_sc[low_sel2]<-spm_low_interp[low_sel2]
        }
        
        if (length(which(flux_sc==0 & specNew$wave>specNew$wave[which(flux_sc==max(flux_sc,na.rm=T))])>0)) {
            highCut<-min(specNew$wave[which(flux_sc==0 & specNew$wave>specNew$wave[which(flux_sc==max(flux_sc,na.rm=T))])])
            high_sel<-which(specNew$wave>=(highCut-500))
            spm_high_sel<-which(spm_wave>=(highCut-500))
            spm_high_interp<-approx(spm_wave,spm_flux_sc,specNew$wave[high_sel])$y
            off_high<-median(spm_high_interp[which(flux_sc[high_sel]>0)]-flux_sc[high_sel[which(flux_sc[high_sel]>0)]],na.rm=T)
            spm_high_interp<-spm_high_interp-off_high
            high_sel2<-which(specNew$wave>=(highCut-10))
            flux_sc[high_sel2]<-approx(specNew$wave[high_sel],spm_high_interp,specNew$wave[high_sel2])$y
        }

  
  
        ##########################################

        wavefac=1e-10
        c<-299792458

        filter=getfilt(band)[,2:3]      
        filt_wave <- filter[,1]
        filt_trans <- filter[,2]

        interp<-approx(filt_wave, filt_trans, specNew$wave)
        filt_trans_interp <- interp$y
        filt_trans_interp[is.na(filt_trans_interp)] = 0

        fluxnu=(wavefac*flux_sc*specNew$wave^2)/c

        flux_conv <- fluxnu*filt_trans_interp

        temp <- fluxnu*filt_trans_interp*specNew$wave
        temp2 <- filt_trans_interp*specNew$wave
        flux_filt <- sum(temp[which(is.finite(temp)==T)])/sum(temp2[which(is.finite(temp2)==T)])
        flux_need <- 10.^(-0.4*(mag[i]+48.6))

        flux_sc_Hz <- fluxnu*(flux_need/flux_filt)
        flux_sc <- (flux_sc_Hz*(2.998e18))/(specNew$wave^2)
        
        name<-paste(outName,'/',outName,'_',id[i],'.fits',sep='')
               

        ### Save spectrum ###
        writeFITSim(cbind(flux_sc), file = name, axDat = axDat, type='single')
        write.fitskey('OBJECT',paste(outName,'_',id[i], sep=''), name, comment = 'redshift', hdu = 1)            
        write.fitskey('Z',z[i], name, comment = 'redshift', hdu = 1)
        write.fitskey('BAND',band, name, comment = 'band', hdu = 1)
        write.fitskey('MAG',mag[i], name, comment = 'magnitude in band', hdu = 1)
        write.fitskey('CRTYPE1','WAVE', name, comment = 'X-axis type', hdu = 1)
        write.fitskey('BUNIT','erg/s/cm2/Angstrom', name, comment = 'flux density unit', hdu = 1)
        write.fitskey('ABMAG AB',mag[i], name, comment = 'SDSS r-mag', hdu = 1)
        write.fitskey('RESOLUTN',2, name, comment = 'resolution', hdu = 1)

      
        

    }
	


}
