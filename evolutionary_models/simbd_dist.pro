function simbd_dist, mass=mass, age=age, sigage=sigage,d=d,sigd=sigd,filter=filter,min_eff=min_eff,n_trials=n_trials,save_plot=save_plot,mdusty=mdusty;,acc_eff=acc_efficiency,acc_eff=acc_eff_sigma

;this script generates a mass density distribution given the following inputs
; mass = a single value in terms of Jupiter masses
; age = a distribution (Gaussian typically)
; distance  = a distribution (Gaussian typically)
; ...
; and one more for hot / cold start interpolation...
; accretion_efficiency = a uniformly distributed value from min_eff to 1

; KW 2018

;=================================================================================

iage=age & id=d	;store inputs as different names

smass=mass*0.000954588	;solar mass units

;n_trials = 10000.	;the more trials the better samples the distribution, but also the longer the routine will take. 

;this loop generates a synthetic photometric probability density distribution 

;test_filter='Mh'
;min_eff=0.5 ;minimum efficiency for core accretion

mags=[] ;this array will store the distribution of generated magnitudes
dists=[]
ages=[] ;storing inputs
	Seed1=systime(/Seconds)
for kkk=0, n_trials do begin

	;generate three randomly distributed gaussian numbers
	Seed=(Seed1-systime(/Seconds))*10000.
	random=randomn(Seed,3)
	random[2]=min_eff+( (1.-min_eff) * randomu(Seed)) ; uses a uniform distribution for the accretion efficiency
	
	;scale the random numbers by the unertainties
	random[0]=random[0]*sigage & random[1]=random[1]*sigd 
	
	;create newly realized parameters
	newage=double(iage+random[0]) & newd=double(id+random[1])  & acceff=double(random[2]) ;& mass=imass

	;ensure that ages are gt 1 Myr and lt 10 Gyr
	;newage=max([newage,1])
	;newage=min([newage,10000])

	if newage gt 10000 or newage lt 1 then begin mag=!values.f_nan & newage=!values.f_nan & newd=!values.f_nan & endif
	
	;print, mass, newage, newd
	;print, random[0]

	;generate the apparent magnitude using Daniel's script
	if not keyword_set(mdusty) then mag=simbd(mass=smass,age=newage,filter=filter)+distmod(newd) $
		else mag=simbd(mass=smass,age=newage,filter=filter,/mdusty)+distmod(newd)
	;print, mag
	
	;print, 'Mass = '+string(mass)+ ' Mag = '+string(mag-distmod(newd))+' Age = '+string(newage)+' D='+string(newd)
	;include efficiency

	if mass le 13. then begin
		;go to flux units
		flux=10.^(mag/2.5)
		fluxeff=(acceff)^(4.0) ;acceff ~ dE/dS ~ T
		;go back to magnitudes
		magdiff=-2.5*ALOG10(fluxeff)
		mag=mag+magdiff
	endif

	mags=[mags,mag]	;store the magnitude in the array
	ages=[ages,newage]
	dists=[dists,newd]

endfor


;remove outliers?
;if n_elements(mags ge 10) then begin
;mags=mags[sort(mags)]
;mags=mags[5:n_elements(mags)-1-5]
;endif

;replace all elements with the median
;mags[*]=median(mags)

;print, mags


;if 1 eq 0 then begin

if keyword_set(save_plot) then begin
if total(finite(mags)) gt 0 then begin

!p.multi=[0,2,2]
	cgps_open, '~/Desktop/inputs.ps'
mags=mags[where(finite(mags))]
   cgHistoplot, mags, BinSize=(max(mags)-min(mags))/50., /Fill,/frequency,title='Photometric Probability Density for Mass = '+string(sigfig(mass,2))+' MJup',charsize=0.7;, Output='~/Desktop/newtest.jpg'

	;cgPS_Close,/nofix;,/png
;endif

;if min(finite(mags)) gt 0 then begin
	;cgps_open, '~/Desktop/Distance.ps'
dists=dists[where(finite(dists))]
   cgHistoplot, dists, BinSize=(max(dists)-min(dists))/50., /Fill,/frequency,title='Distance Prior',charsize=0.7;, Output='~/Desktop/newtest.jpg'

	;cgPS_Close,/nofix;,/png
;endif

;if min(finite(mags)) gt 0 then begin
	;cgps_open, '~/Desktop/Age.ps'
ages=ages[where(finite(ages))]
   cgHistoplot, ages, BinSize=(max(ages)-min(ages))/50., /Fill,/frequency,title='Age Prior',charsize=0.7;, Output='~/Desktop/newtest.jpg'

	cgPS_Close,/nofix;,/png

	;hak
endif
endif
;endif
return, mags

;save,filename='~/Desktop/mags.sav',mags


end
