get.spec=function(file, row=1, logw=F, sc=1, xunit='ang', yunit='ang', z=NA, RA=NA, DEC=NA) { 
  specf<-readFITS(file)  
  if (length(as.numeric(specf$hdr[which(specf$hdr=='Z')+1])) >0 & is.finite(z)==F) {z<-as.numeric(specf$hdr[which(specf$hdr=='Z')+1])}
  if (length(as.numeric(specf$hdr[which(specf$hdr=='RA')+1])) >0 & is.finite(RA)==F) {RA<-as.numeric(specf$hdr[which(specf$hdr=='RA')+1])}
  if (length(as.numeric(specf$hdr[which(specf$hdr=='DEC')+1])) >0 & is.finite(DEC)==F) {DEC<-as.numeric(specf$hdr[which(specf$hdr=='DEC')+1])}
  if (length(dim(specf$imDat))>1) {flux <- specf$imDat[,row]*sc}
  if (length(dim(specf$imDat))==1) {flux <- specf$imDat*sc}
  CRVAL<-as.numeric(specf$hdr[which(specf$hdr=='CRVAL1')+1])
  CRPIX<-as.numeric(specf$hdr[which(specf$hdr=='CRPIX1')+1]) 
  CDELT<-as.numeric(specf$hdr[which(specf$hdr=='CDELT1')+1])
  if (length(CDELT)==0) {CDELT<-as.numeric(specf$hdr[which(specf$hdr=='CD1_1')+1])}
  
  wave <- (CRVAL+(c(0:(specf$axDat$len[1]-1))*CDELT))-(CDELT*(CRPIX-1))
  if (logw==T) {
    wave<-wave-CDELT
    wave<-10.^wave
  }  
  spec <- list(wave,flux,xunit,yunit, z, RA, DEC)
  names(spec) <- c('wave', 'flux', 'xunit', 'yunit', 'z','RA','DEC')
  return=spec      
}
