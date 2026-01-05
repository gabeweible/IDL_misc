function bd_massdd,objmag=objmag,sigmag=sigmag,filter=filter,age=age,sigage=sigage,$
d=d,sigd=sigd,min_eff=min_eff,minmass=minmass,maxmass=maxmass,stepmass=stepmass,$
n_trials=n_trials,mdusty=mdusty;,acc_eff=acc_efficiency,acc_eff=acc_eff_sigma

;first step: generate a photometric probability density distribution for a range of masses
if not keyword_set(minmass) then minmass=0.
if not keyword_set(maxmass) then maxmass=80. ;Jupiter masses
if not keyword_set(stepmass) then stepmass=1. 
if not keyword_set(n_trials) then n_trials=1000
probs=[] & masses=[];where to store the probabilities for each mass bin

n_masses = FIX((maxmass-minmass) / stepmass)

for mass_i = 0,n_masses do begin

	mass = minmass + mass_i*stepmass

;for mass=minmass,maxmass,stepmass do begin	
	;print, 'Testing mass ', mass
	if not keyword_set(mdusty) then mags=simbd_dist(mass=mass,age=age,sigage=sigage,d=d,sigd=sigd,min_eff=min_eff,n_trials=n_trials,filter=filter)
	if keyword_set(mdusty) then mags=simbd_dist(mass=mass,age=age,sigage=sigage,d=d,sigd=sigd,min_eff=min_eff,n_trials=n_trials,filter=filter,/mdusty);,
	

	;what is the probability that the object is this mass?
	;minmag=objmag-(1.*sigmag) & maxmag=objmag+(1.*sigmag)
	;prob=double(n_elements(where(mags ge minmag and mags le maxmag)));/float(n_elements(mags))
	;if n_elements(where(mags ge minmag and mags le maxmag)) eq 1 then if where(mags ge minmag and mags le maxmag) eq -1 then prob=0
	

;find distances
	;dists=mags-objmag

	;make a guassian
	parms=[1.,objmag,sigmag] & gauss=gaussian(mags,parms)

	;calculate relative probability (arbitrary normalization) as the sum of the values of the gaussian function at the simualted magnitudes
	prob=total(gauss)


	;if n_elements(where(mags ge minmag and mags le maxmag)) eq 1 then if where(mags ge minmag and mags le maxmag) eq -1 then prob=0
	
	;prior=1.	;set the mass prior here, uniform = 1
	
	;prob=prob*prior
	probs=[probs,prob]
	masses=[masses,mass]
endfor

	;ensure total probability is unity

	;print, minmag, maxmag
	;print, masses
	;print, probs
	;plot, masses, probs
	;hak
	
	probs=probs/total(probs)


;now the probability that the object is each mass is given in the prob array


	
return,probs
end
