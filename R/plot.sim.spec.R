plot.sim.spec=function(simName, id){
  for (i in 1:length(id)){
    spec<-get.spec(paste(simName,'/',simName,'_',id[i],'.fits', sep=''))
    plot.spec(spec, main=paste(simName,'_',id[i],'.fits', sep=''))
    if (length(id)>1){readline('Plot Next Spectrum?')}
  }
}
