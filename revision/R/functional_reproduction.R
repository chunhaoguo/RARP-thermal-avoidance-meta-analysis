options(stringsAsFactors=FALSE)
root <- normalizePath('.',mustWork=TRUE)
dir.create('results',showWarnings=FALSE)
d <- read.delim('data/inputs/functional_trajectory_source.tsv',na.strings=c('', 'NA'))
d$effect <- d$low <- d$high <- d$variance <- NA_real_
for (i in seq_len(nrow(d))) {
  if (d$panel[i]=='iief5') {
    est <- d$athermal_mean[i]-d$energy_mean[i]
    v <- d$athermal_sd[i]^2/d$athermal_n[i]+d$energy_sd[i]^2/d$energy_n[i]
    ci <- est+c(-1,1)*qnorm(.975)*sqrt(v)
  } else {
    est <- (d$athermal_events[i]/d$athermal_total[i])/(d$energy_events[i]/d$energy_total[i])
    v <- 1/d$athermal_events[i]-1/d$athermal_total[i]+1/d$energy_events[i]-1/d$energy_total[i]
    ci <- exp(log(est)+c(-1,1)*qnorm(.975)*sqrt(v))
  }
  d$effect[i]<-est;d$variance[i]<-v;d$low[i]<-ci[1];d$high[i]<-ci[2]
}
expected <- read.delim('data/inputs/expected_functional_effects.tsv')
stopifnot(identical(d$panel,expected$panel),identical(d$time_label,expected$time_label))
for (v in c('effect','low','high','variance')) stopifnot(max(abs(d[[v]]-expected[[v]]))<1e-10)
d$analysis_role <- ifelse(d$panel=='iief5','exploratory conditional reconstruction','individual study effect')
d$time_origin <- ifelse(d$study_label=='Fuschi 2025','catheter removal','study-specific postoperative follow-up')
write.table(d,'results/functional_effects_verified.tsv',sep='\t',row.names=FALSE,quote=FALSE)

# Protocol-triggered two-study non-randomised blood-loss calculation. Exploratory only;
# the randomised Fuschi estimate is displayed separately under the design-stratified rule.
p <- read.delim('results/perioperative_effects.tsv')
p <- p[p$outcome=='estimated_blood_loss' & p$design=='Non-randomized',]
y<-p$effect;v<-p$variance
obj<-function(t){w<-1/(v+t);mu<-sum(w*y)/sum(w);.5*(sum(log(v+t))+log(sum(w))+sum(w*(y-mu)^2))}
tau<-optimize(obj,c(0,10000))$minimum
wr<-1/(v+tau);mu<-sum(wr*y)/sum(wr);se<-sqrt(1/sum(wr))
wf<-1/v;mf<-sum(wf*y)/sum(wf);sf<-sqrt(1/sum(wf));Q<-sum(wf*(y-mf)^2)
b<-data.frame(model=c('REML Wald exploratory','Common effect exploratory'),MD=c(mu,mf),low=c(mu-1.96*se,mf-1.96*sf),high=c(mu+1.96*se,mf+1.96*sf),tau2=c(tau,0),I2=100*max(0,(Q-1)/Q))
write.table(b,'results/blood_loss_exploratory.tsv',sep='\t',row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),'results/sessionInfo.txt')
cat('All 11 functional effects match the archived estimates within 1e-10.\n')
