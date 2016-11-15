plot.spec=function(spec, z=NA, plot_lines=T, plot_AGN_lines=F, oplot=F, col='black', xlim=NA, ylim=NA, main=NA) {
  
  xr<-c(min(spec$wave),max(spec$wave))
  yr<-c(min(spec$flux[is.finite(spec$flux)]),max(spec$flux[is.finite(spec$flux)]))
  
  if (is.finite(xlim[1])==T) {xr=xlim}
  if (is.finite(ylim[1])==T) {yr=ylim}
  
  if (is.finite(z)==F) {z<-spec$z[1]}
  if (spec$xunit[1]=='ang') {xtit<-'Wavelength, Ang'}
  if (spec$xunit[1]=='hz') {xtit<-'Frequency, Hz'}
  if (spec$xunit[1]=='m') {xtit<-'Wavelength, m'}
  if (spec$xunit[1]=='micron') {xtit<-'Wavelength, micron'}
  if (spec$xunit[1]=='nm') {xtit<-'Wavelength, nm'}
  if (spec$yunit[1]=='ang') {ytit<-'Flux, ergs/sec/cm^2/ang'}
  if (spec$yunit[1]=='hz') {ytit<-'Flux, ergs/sec/cm^2/Hz'}
  if (spec$yunit[1]=='Jy') {ytit<-'Flux, Jy'}
  
  lines<-load.lines()   
  line_names <- as.character(lines$name)
  
  lines_AGN<-load.lines.AGN()
  line_names_AGN<- as.character(lines_AGN$name)
  
  if (spec$xunit[1]=='ang') {line_x <- as.numeric(lines$wave_ang)*(1+z)}
  if (spec$xunit[1]=='hz') {line_x <- as.numeric(lines$freq_hz)/(1+z)}
  if (spec$xunit[1]=='micron') {line_x <- as.numeric(lines$wave_micron)*(1+z)}
  if (spec$xunit[1]=='m') {line_x <- as.numeric(lines$wave_m)*(1+z)}
  if (spec$xunit[1]=='nm') {line_x <- as.numeric(lines$wave_nm)*(1+z)}
  if (spec$xunit[1]=='ang') {line_x_AGN <- as.numeric(lines_AGN$wave_ang)*(1+z)}
  if (spec$xunit[1]=='hz') {line_x_AGN <- as.numeric(lines_AGN$freq_hz)/(1+z)}
  if (spec$xunit[1]=='micron') {line_x_AGN <- as.numeric(lines_AGN$wave_micron)*(1+z)}
  if (spec$xunit[1]=='m') {line_x_AGN <- as.numeric(lines_AGN$wave_m)*(1+z)}
  if (spec$xunit[1]=='nm') {line_x_AGN <- as.numeric(lines_AGN$wave_nm)*(1+z)}
  
  
  if (oplot==F) {
    if (length(dev.list())==0) {quartz(width=12, height=6)}
    magplot(spec$wave, spec$flux, type='l', xlab=xtit, ylab=ytit, col=col, xlim=xr, ylim=yr, main=main)
    if (plot_lines==T) {
      for (i in 1:length(line_names)) {
        if (lines$stellar[i]==F) {
          lines(c(line_x[i], line_x[i]),c(-1e-100,1e100), col='blue', lty=2)
          text(line_x[i], max(spec$flux), line_names[i], col='blue',cex=0.4)
        }
        if (lines$stellar[i]==T) {
          lines(c(line_x[i], line_x[i]),c(-1e-100,1e100), col='darkgreen', lty=2)
          text(line_x[i], max(spec$flux), line_names[i], col='darkgreen',cex=0.4)
        }
      }
    }
    if (plot_AGN_lines==T) {
      for (i in 1:length(line_names_AGN)) {
        lines(c(line_x_AGN[i], line_x_AGN[i]),c(-1e-100,1e100), col='orange', lty=2)
        text(line_x_AGN[i], max(spec$flux), line_names_AGN[i], col='orange',cex=0.4)           
      }
    }
    
  }
  if (oplot==T) {lines(spec$wave, spec$flux, col=col)}
}
