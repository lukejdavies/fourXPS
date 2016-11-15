load.lines.AGN=function(){
  c<-299792458
  names <- c('NV','CIV', 'CIII', 'FeII', 'MgII' )
  stellar <- c(F,F,F,F,F)
  wave_ang <- c(1240,1549,1909, 2500, 2798)
  freq_hz <- c/(wave_ang/(10.^10))
  wave_m <- wave_ang/(10.^10)
  wave_micron <- wave_ang/(10.^4)
  wave_nm <- wave_ang/(10.^1)
  lines<-data.frame(names, wave_ang,wave_m,wave_micron,wave_nm, freq_hz, stellar)
  return(lines) 
}