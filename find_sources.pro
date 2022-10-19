;This routine automatically finds sources above a particular SNR sigma in an ADI-type dataset. i.e. applicable when the noise profile is more uniform when accounted for at a particular radius away from the center than among the frame as a whole. 



;   0         None           No scaling whatsoever is done.
;   1         Linear         scaled = BytScl(image, MIN=minValue, MAX=maxValue)
;   2         Clip           A histogram stretch, with a percentage of pixels clipped at both the top and bottom
;   3         Gamma          scaled = GmaScl(image, MIN=minValue, MAX=maxValue, Gamma=gamma)
;   4         Log            scaled = LogScl(image, MIN=minValue, MAX=maxValue, Mean=mean, Exponent=exponent)
;   5         Asinh          scaled = AsinhScl(image, MIN=minValue, MAX=maxValue, Beta=beta)
;   6         SquareRoot     A linear stretch of the square root histogram of the image values.
;   7         Equalization   A linear stretch of the histogram equalized image histogram.
;   8         Gaussian       A Gaussian normal function is applied to the image histogram.
;   9         MODIS          Scaling done in the differential manner of the MODIS Rapid Response Team
;                            and implemented in the Coyote Library routine ScaleModis.
;   10        StdDev         Standard deviation stretch. scaled = SDevScl(image, MULTIPLIER=2).

;Kevin Wagner - Steward Observatory

pro find_sources,frame,reference=reference,sigma=sigma,correction_factor=correction_factor $
,FWHM=FWHM,clip_level=clip_level,platescale=platescale,inrad=inrad,outrad=outrad,stretch=stretch,$
minvalue=minvalue,maxvalue=maxvalue,no_2d=no_2d,no_smooth=no_smooth,ct=ct,filter=filter


;tunable parameters

if not keyword_set(correction_factor) then correction_factor=1.0
if not keyword_set(sigma) then sigma=5. 
;if not keyword_set(splatescale) then platescale=1.
if not keyword_set(FWHM) then FWHM=4.;(wavelength)/(diameter) * 206265. / platescale
aper_rad=FWHM/2.
if not keyword_set(clip_level) then clip_level=0.5

if not keyword_set(reference) then do_phot=0 else do_phot=1
if not keyword_set(stretch) then stretch='linear'

if keyword_set(no_2d) then do_2d=0 else do_2d=1	;2D SNR and contrast maps take a while to generate

if do_phot then begin
	reference=readfits(reference)

	rsizes=size(reference) & rxsize=rsizes[1] & rysize=rsizes[2]
	
	aper, reference,rxsize/2.,rysize/2., ref_flux,fluxerr,sky,skyerr,1.75,aper_rad,sky_rad,[-99E99,99E99],/silent,/flux,SETSKYVAL=0,/exact


endif; do_phot if
if not do_phot then ref_flux=1
;end user input, begin script

;print, 'Beginning tests with aperture radius = (px)',aper_rad




sky_rad=[5,7] ;does not do anything, just needs to be set

;if not keyword_set(frame) then frame='/home/kevin/Desktop/LBTI_MWC758_Lp_combined.fits' ;to be deleted


pos=strpos(frame,'/',/reverse_search)
location=strmid(frame,0,pos+1)
name=strmid(frame,pos+1,strlen(frame)-1)

frame=readfits(frame)

frame[where(frame le 0. and frame ge 0.)]=!values.f_nan



;close,1; commented out because of writing to .txt error
;openw,1,location+name+ '_sources.txt',width=400; commented out because of writing to .txt error





;read in an frame with the input "frame" paramater (expects a 2d array)

;find x,y dimensions (xsize,ysize)

sizes=size(frame) & xsize=sizes[1] & ysize=sizes[2]

;make smaller
;print, outrad
if not keyword_set(outrad) then outrad=xsize/2.
;print, xsize/2.-outrad-16.,xsize/2.+outrad-1+16.
if outrad lt xsize/2.-16. then begin
	frame=frame[xsize/2.-outrad-16.:xsize/2.+outrad-1+16.,ysize/2.-outrad-16.:ysize/2.+outrad-1+16.]
	sizes=size(frame) & xsize=sizes[1] & ysize=sizes[2]
endif


;loop through radii from 1 to rsize
rsize=min([xsize,ysize])/2.;sqrt( (xsize/2.)^2. + (ysize/2.)^2 )

if not keyword_set(inrad) then inrad=1.
if not keyword_set(outrad) then outrad=(xsize/2.) - 4.*fix(FWHM)

;print, 'frame has radius =',rsize


	;window,0,xs=xsize,ys=ysize
	;w=window(dimensions=[xsize,ysize])
	loadct,36,/silent

		frame[where(finite(frame) eq 0)]=0.

		disp_frame=frame

		;disp_frame[where(disp_frame gt max(disp_frame)/disp_contrast)]=max(disp_frame)/disp_contrast

		;disp_frame= disp_frame < maxval
		;disp_frame=disp_frame+min(disp_frame)
		
		;if disp_type eq 'linear' then tv,disp_frame
		;if disp_type eq 'sqrt' then tvscl,sqrt(disp_frame)
		;if disp_type eq 'log' then tvscl,alog(disp_frame)
		;if disp_type eq 'asin' then tvscl,asin(disp_frame)
		;if disp_type eq 'sin' then tvscl,sin(disp_frame)

		;if disp_type eq 'squared' then tvscl,-(disp_frame)*(disp_frame)
		
		;graphic=image(disp_frame,/no_toolbar,image_dimensions=[xsize,ysize],margin=0,dimensions=[600,600],RGB_TABLE=37)

;mydevice = !D.NAME
; Set plotting to PostScript:
;SET_PLOT, 'PS'
;DEVICE, FILENAME=location+'_source_map.ps',/PORTRAIT, /COLOR,/INCHES, YSIZE=3.

disp_frame=frame

for xx=0,xsize-1 do for yy=0,ysize-1 do if sqrt((xx-(xsize/2.))^2. + (yy-(ysize/2.))^2.) lt inrad or $
	sqrt((xx-(xsize/2.))^2. + (yy-(ysize/2.))^2.) gt outrad then disp_frame[xx,yy]=!values.f_nan

disp_frame[where(disp_frame ge 0 and disp_frame le 0)]=!values.f_nan



sources=[]
sources_x=[0.]
sources_y=[0.]
sources_sigma=[]
sources_contrast=[]
xlast=0 & ylast=0 & flipswitch=0 & stopswitch=0

for rr=inrad,outrad do begin

	;start at north, progress by each pixel (to 2PiR) and calculate flux in each aperture
	ntheta=(2.*!PI*float(rr))/(fix(FWHM)) ;about one angular bin per resolution element

	;print, '--- Searching for sources at radius (px): ',rr
	;print, '--- Apertures at this radius: ',ntheta

	for itheta=0.,ntheta do begin
		theta=(float(itheta)/float(ntheta))*2.*!PI;+(!PI/31.)
		xtest=float(xsize)/2.+(float(rr)*cos(theta)) & ytest=float(ysize)/2.+(float(rr)*sin(theta))
		;print, size(frame)
		aper, frame, xtest, ytest, flux,fluxerr,sky,skyerr,1.75,aper_rad,sky_rad,[-99E99,99E99],/silent,/flux,SETSKYVAL=0,/exact

		if itheta eq 0 then fluxarr=flux else if flux ne 0 then fluxarr=[fluxarr,flux]
		if itheta eq 0 then xarr=xtest else xarr=[xarr,xtest]
		if itheta eq 0 then yarr=ytest else yarr=[yarr,ytest]
	endfor ;;itheta for

	;print, fluxarr
	;print, 'Standard deviation of flux measure: ',stddev(fluxarr,/NAN)
	;print, 'Mean of abs of flux measures: ',mean(abs(fluxarr),/NAN)





	tfluxarr=fluxarr[where(finite(fluxarr) eq 1)]

	;tfluxarr=fluxarr[where(fluxarr le 0)] ;take the negatives so that real sources do not confuse the selection.
	;tfluxarr=fluxarr ;take all of it

	txarr=xarr;[where(fluxarr le 0)]
	tyarr=yarr;[where(fluxarr le 0)]

	fluxarrl=fluxarr[where(txarr le xsize/2.)]
	maxtestl=max(abs(fluxarrl),/NAN)	;now re-computes the noise
	;tfluxarrl=fluxarrl[where(fluxarrl le 0)]

	tfluxarrl=fluxarrl[where(finite(fluxarrl) eq 1)]

	fluxarrr=fluxarr[where(txarr gt xsize/2.)]
	maxtestr=max(abs(fluxarrr),/NAN)	;now re-computes the noise
	;tfluxarrr=fluxarrr[where(fluxarrr le 0)]

	tfluxarrr=fluxarrr[where(finite(fluxarrr) eq 1)]




	fluxarrt=fluxarr[where(tyarr le xsize/2.)]
	maxtestt=max(abs(fluxarrt),/NAN)	;now re-computes the noise
	;tfluxarrl=fluxarrl[where(fluxarrl le 0)]

	tfluxarrt=fluxarrt[where(finite(fluxarrt) eq 1)]

	fluxarrb=fluxarr[where(tyarr gt xsize/2.)]
	maxtestb=max(abs(fluxarrb),/NAN)	;now re-computes the noise
	;tfluxarrr=fluxarrr[where(fluxarrr le 0)]

	tfluxarrb=fluxarrb[where(finite(fluxarrb) eq 1)]



	test=stddev(tfluxarr,/NAN)	;if taking negatives multiply by two because now we have taken half of the distribution (assuming it is distributed about zero)
	;test=mean(abs(tfluxarr),/NAN)	

testl=stddev(tfluxarrl,/NAN)
testr=stddev(tfluxarrr,/NAN)


testt=stddev(tfluxarrt,/NAN)
testb=stddev(tfluxarrb,/NAN)


;maxtestl=max(abs(fluxarrl),/NAN)
;maxtestr=max(abs(fluxarrr),/NAN)



;if testl/testr gt 4 then test=testr
;if testr/testl gt 4 then test=testl


	tfluxarr=fluxarr[where(abs(fluxarr) lt 5.*test)] ;focus on just the noise elements, ignores real signals and theit negative counterparts from post processing
tfluxarrl=fluxarrl[where(abs(fluxarrl) lt 5.*testl)]
tfluxarrr=fluxarrr[where(abs(fluxarrr) lt 5.*testr)]



tfluxarrt=fluxarrl[where(abs(fluxarrt) lt 5.*testt)]
tfluxarrb=fluxarrr[where(abs(fluxarrb) lt 5.*testb)]



	test=stddev(tfluxarr,/NAN)	;now re-computes the noise
	testmean=mean(tfluxarr,/NAN)
	test=test*sqrt(1.+(1./n_elements(tfluxarr)))	;correction for small sample statistics

	;adding a protocol for bright binaries.
	testl=stddev(tfluxarrl,/NAN)	;now re-computes the noise
	testmeanl=mean(tfluxarrl,/NAN)
	testl=testl*sqrt(1.+(1./n_elements(tfluxarrl)))
	testr=stddev(tfluxarrr,/NAN)	;now re-computes the noise
	testmeanr=mean(tfluxarrr,/NAN)
	testr=testr*sqrt(1.+(1./n_elements(tfluxarrr)))


	testt=stddev(tfluxarrt,/NAN)	;now re-computes the noise
	testmeant=mean(tfluxarrt,/NAN)
	testt=testt*sqrt(1.+(1./n_elements(tfluxarrt)))
	testb=stddev(tfluxarrb,/NAN)	;now re-computes the noise
	testmeanb=mean(tfluxarrb,/NAN)
	testb=testb*sqrt(1.+(1./n_elements(tfluxarrb)))

	testn=n_elements(tfluxarr)
	testnl=n_elements(tfluxarrl)
	testnr=n_elements(tfluxarrr)

	testnt=n_elements(tfluxarrt)
	testnb=n_elements(tfluxarrb)
if testl/testr gt 3. and n_elements(tfluxarrr) gt 4 then test=testr
if testr/testl gt 3. and  n_elements(tfluxarrl) gt 4 then test=testl


if testl/testr gt 3. and n_elements(tfluxarrr) gt 4 then testn=testnr
if testr/testl gt 3. and  n_elements(tfluxarrl) gt 4 then testn=testnl


if testl/testr gt 3. and n_elements(tfluxarrr) gt 4 then testmean=testmeanr
if testr/testl gt 3. and  n_elements(tfluxarrl) gt 4 then testmean=testmeanl


if testt/testb gt 3. and n_elements(tfluxarrb) gt 4 then test=testb
if testb/testt gt 3. and  n_elements(tfluxarrt) gt 4 then test=testt


if testt/testb gt 3. and n_elements(tfluxarrb) gt 4 then testn=testnb
if testb/testt gt 3. and  n_elements(tfluxarrt) gt 4 then testn=testnt


if testb/testt gt 3. and n_elements(tfluxarrb) gt 4 then testmean=testmeant
if testt/testb gt 3. and  n_elements(tfluxarrt) gt 4 then testmean=testmeanb

	;test=mean(abs(tfluxarr),/NAN)

	;print, where(tfluxarr ge 0 and tfluxarr le 0)
	;if rr gt inrad + 5 then if n_elements(where(tfluxarr ge 0 and tfluxarr le 0)) gt 1 then test=mean([testlast,testlastsq,testlastcu,testlastfo,testlastfi])	;helps for corners

	;if rr gt inrad+4 then testlastfi=testlastfo
	;if rr gt inrad+3 then testlastfo=testlastcu
	;if rr gt inrad+2 then testlastcu=testlastsq
	;if rr gt inrad+1 then testlastsq=testlast
	;testlast=test

;helps with corners

;posframe=frame > 0

;print, rr-inrad
;print, (outrad-inrad/2.)
;print, outrad,inrad
if rr-inrad gt (outrad-inrad)/2. then begin
	fit=poly_fit(rs[(outrad-inrad)/2.:rr-inrad-1],tests[(outrad-inrad)/2.:rr-inrad-1],0)
	if n_elements(where(tfluxarr ge 0 and tfluxarr le 0)) gt 1 then $
	 test=fit[0]
endif


if rr eq inrad then testns=testn else testns=[testns,testn]
if rr eq inrad then tests=test else tests=[tests,test]
if rr eq inrad then testmeans=testmean else testmeans=[testmeans,testmean]
if rr eq inrad then rs=rr else rs=[rs,rr]

	for ii=0,n_elements(fluxarr)-1 do begin
		if (fluxarr[ii]-testmean)/test ge sigma/2. then begin

				;smframe=frame > 0.
			xcntr=xarr[ii] & ycntr=yarr[ii]
			cntrd,frame,xarr[ii],yarr[ii],xcntr,ycntr,FWHM,/silent;,extendbox=3.*FWHM

			;if xcntr eq -1 or ycntr eq -1 then begin

				;smframe=smooth(frame,3.)
			;	cntrd,frame,xarr[ii],yarr[ii],xcntr,ycntr,FWHM,extendbox=FWHM;	,/debug

			;endif

			;if xcntr eq -1 or ycntr eq -1 then begin




			;endif


			
			if xcntr gt -1 and ycntr gt -1 then begin
			;if xcntr eq -1  then xcntr=xarr[ii]
			;if ycntr eq -1  then ycntr=yarr[ii]
			;repeat measurement
			aper, frame, xcntr,ycntr, flux,fluxerr,sky,skyerr,1.75,aper_rad,sky_rad,[-99E99,99E99],/silent,/flux,SETSKYVAL=0,/exact

			if (flux-testmean)/test gt sigma then begin


			;if xcntr lt 0 then xcntr=xarr[ii]
			;if ycntr lt 0 then ycntr=yarr[ii]
	

			;test local significance
			padsz=24
			;padsz=2
			boxsz=min([padsz,4.*fix(FWHM)])
			pad=max([padsz,boxsz])
			
			if xcntr-pad gt 0 and xcntr+pad-1. lt xsize-1 and ycntr-pad gt 0 and ycntr+pad-1 lt ysize-1. then begin
	
			box=frame[xcntr-boxsz:xcntr+boxsz-1,ycntr-boxsz:ycntr+boxsz-1]
			box[max([0,boxsz/2.-(1.5*fix(FWHM))]):min([2.0*xcntr,boxsz/2.+(1.5*fix(FWHM))]),max([0,boxsz/2.-(1.5*fix(FWHM))]):min([2.0*ycntr,boxsz/2.+(1.5*fix(FWHM))])]=!values.f_nan
			box[where(box eq 0.)]=!values.f_nan
			;box=box[where(box gt 0)]

			boxstdv=stddev(box,/NAN)

			aaa=max(frame[xcntr-FWHM:xcntr+FWHM,ycntr-FWHM:ycntr+FWHM]) - mean(box,/NAN)
			endif else begin
				box=0.
				boxstdv=1000.

				aaa=0 ;should essentially ignore edges with this and the next if statement
			endelse
			if rr lt 10.*FWHM then boxstdv=0;exclude 10 l/d since local stdv is varying rapidly within here 

			if aaa gt min([sigma,2])*boxstdv then begin	

			if sources_x eq !null then flipswitch=0
			if sources_x eq !null then sources_x=0.
			if sources_y eq !null then sources_y=0.

			x_dist=sources_x-xcntr
			y_dist=sources_y-ycntr
			dists=sqrt(  ((x_dist)^2.)  + ((y_dist)^2.)     )
			min_dist= min(dists)

			if  min_dist ge FWHM   then begin	


			

			;print, '----- Centered source significance = ',(flux-testmean)/test
			if do_phot then print, '----- Source contrast = ',flux/ref_flux
		
				if flipswitch ne 1 then begin sources_x=[] & sources_y=[] & endif & flipswitch=1

							smframe=frame > 0.
			
			
			;print, '--- Found a significant source!'
			;print, '----- Source test significance = ',(fluxarr[ii]-testmean)/test
			;print, '----- Source test x,y = ',xarr[ii],yarr[ii]

			;print, '----- Source centroid x,y = ',xcntr,ycntr

				

			

	;if xcntr-19. ge 0 and xcntr+20. le xsize-1 and ycntr-19. ge 0 and ycntr+20. le ysize-1. then begin
	;print, xcntr, ycntr			
	source=frame[xcntr-19.:xcntr+20.,ycntr-19.:ycntr+20.]

		fit=gauss2dfit(source,A,/tilt)

		;sourcewidth=(A[4]+A[5])/2.
		sourcewidth=min([ a[4],a[5] ])	
		if max([a[4],a[5]]) / min([a[4],a[5]]) lt 1.2 then begin ;test ratio of widths
		if sourcewidth ge 3.5 then begin


			sources=[ [[sources]],[[source]] ]

		;add here a FWHM test
		;fita=[1,1,1,1,0,0,1]
		;A=[0.,0.,0.,0.,sources_x,sources_y,0.]
		
	;endif
			sources_x=[sources_x,xcntr]
			sources_y=[sources_y,ycntr]
			sources_sigma=[sources_sigma,(flux-testmean)/test]



			sources_contrast=[sources_contrast,flux/double(ref_flux)]
			xlast=xcntr & ylast=ycntr



			endif ;local if

		endif ;width if
		endif ;width ratio if
			endif ;min if
			endif ;flux/test if
			endif ;-1 if

			;hak
		endif
	endfor ;ii for

	
	;call the above value the noise
;aper,ifsfluxcube[*,*,i],ifsX,ifsY,flux,errap,sky,skyerr, phpadu, aperture, skyradst,$
	;  stbadpix, /silent, /flux 

endfor ;rr for


;2D sensitivity map
if do_2d then begin



result=frame
result[*]=!values.f_nan
varresult=result
if do_phot then contrast_result=result
;print, 'Generating 2D sensitivity map...'
for xx=float(fix(FWHM)),xsize-fix(FWHM) do begin
for yy=float(fix(FWHM)),ysize-fix(FWHM) do begin
	
	rad=sqrt( (xx-(xsize/2.))^2. + (yy-(ysize/2.))^2. ) 

	if rad gt inrad and rad lt outrad then begin
	testrad=tests[where(round(rs) eq round(rad))]
	testnrad=testns[where(round(rs) eq round(rad))]
	testmeanrad=testmeans[where(round(rs) eq round(rad))]
	aper, frame,xx,yy, flux,fluxerr,sky,skyerr,1.75,aper_rad,sky_rad,[-99E99,99E99],/silent,/flux,SETSKYVAL=0,/exact

	result[xx,yy]=(flux-mean(testmeanrad))/(testrad)
	varresult[xx,yy]=(flux-mean(testmeanrad))/(testrad*testrad)

	if do_phot then contrast_result[xx,yy]=correction_factor*flux/ref_flux
	endif

endfor
endfor

varresult[where(frame le 0. and frame ge 0.)]=!values.f_nan
if not keyword_set(no_smooth) then varresult=smooth(varresult,FWHM,/NAN)

result[where(frame le 0. and frame ge 0.)]=!values.f_nan
if not keyword_set(no_smooth) then result=smooth(result,FWHM,/NAN)
if do_phot then contrast_result[where(frame le 0. and frame ge 0.)]=!values.f_nan

writefits,location+name+ '_SNR_map.fits',result
writefits,location+name+ '_SVR_map.fits',varresult
if do_phot then writefits,location+name +  '_contrast_map.fits',contrast_result

if sources_x ne !null then begin
sources=[]
for jjj=0,n_elements(sources_x)-1 do begin
	if sources_x[jjj]-19. ge 0 and sources_x[jjj]+20. le xsize - 1 and sources_y[jjj]-19. ge 0 and sources_y[jjj]+20. le ysize - 1 then begin			
			source=result[sources_x[jjj]-19.:sources_x[jjj]+20.,sources_y[jjj]-19.:sources_y[jjj]+20.]

			sources=[ [[sources]],[[source]] ]
	endif
endfor

endif

endif




;DEVICE, /CLOSE
;SET_PLOT, mydevice

;minsrcx=sources_x[where(sources_contrast eq min(sources_contrast))]
;minsrcy=sources_y[where(sources_contrast eq min(sources_contrast))]
test_disp_frame=disp_frame
;block inner half of image to get an idea of the image noise without the central star
for xx=0,xsize-1 do for yy=0,ysize-1 do if sqrt( (xx-(xsize/2.))^2. + (yy-(ysize/2.))^2. ) lt xsize/4 then test_disp_frame[xx,yy]=!values.f_nan

disp_framel=disp_frame[0:xsize/2,*]
disp_framer=disp_frame[xsize/2:xsize-1,*]


disp_framet=disp_frame[*,0:ysize/2]
disp_frameb=disp_frame[*,ysize/2:ysize-1]

test_disp_frame[where(finite(test_disp_frame) ne 1)]=0.
test_disp_framel=test_disp_frame[0:xsize/2,*]
test_disp_framer=test_disp_frame[xsize/2:xsize-1,*]

test_disp_framet=test_disp_frame[*,0:ysize/2]
test_disp_frameb=test_disp_frame[*,ysize/2:ysize-1]

sig=abs(stddev(disp_frame[where(test_disp_frame lt 0)],/NAN))
sigl=abs(stddev(disp_framel[where(test_disp_framel lt 0)],/NAN))
sigr=abs(stddev(disp_framer[where(test_disp_framer lt 0)],/NAN))


sigt=abs(stddev(disp_framet[where(test_disp_framet lt 0)],/NAN))
sigb=abs(stddev(disp_frameb[where(test_disp_frameb lt 0)],/NAN))

maxl=max(test_disp_framel)
maxr=max(test_disp_framer)


maxt=max(test_disp_framet)
maxb=max(test_disp_frameb)

sig=stddev(test_disp_frame,/NAN)
sigl=stddev(test_disp_framel,/NAN)
sigr=stddev(test_disp_framer,/NAN)


sigt=stddev(test_disp_framet,/NAN)
sigb=stddev(test_disp_frameb,/NAN)
;test for binaries

if maxl/maxr gt 3. then sig=sigr
if maxr/maxl gt 3. then sig=sigl


if maxt/maxb gt 3. then sig=sigb
if maxb/maxt gt 3. then sig=sigt

sig=sig

minvalue=-5.0*sig
maxvalue=15.0*sig
if do_phot then maxcon=maxvalue/max(reference)

if do_phot then mincon=minvalue/max(reference)
;print, 'Max, min, mean = ',maxvalue,minvalue,meanvalue
;hak

plotcolor='Snow'

;if do_phot then disp_frame=disp_frame/ref_flux

cgps_open,location+name+ '_original.ps'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,xsize,ysize,/device

 ;  image = cgdemodata(18) cgImage,disp_frame,stretch=stretch,clip=clip_level,xrange=[0,xsize],yrange=[0,ysize],/keep_aspect_ratio,minvalue=minvalue,maxvalue=maxvalue,mean=meanvalue;,position=[0,0,1,1]


;begin annotations

if do_phot then range=[0,maxcon] else range=[0,maxvalue]
if do_phot then cbtitle='Contrast' else cbtitle='ADU/Beam'
if do_phot then tickinterval=maxcon/5. else tickinterval=(maxvalue-minvalue)/5.

if total(range) eq 0 then range=[0,max(disp_frame,/nan)]

cgImage,disp_frame,stretch=stretch,clip=clip_level,xrange=[0,xsize],yrange=[0,ysize],/keep_aspect_ratio,minvalue=minvalue,maxvalue=maxvalue,mean=meanvalue;,position=[0,0,1,1]

;print, range
!p.psym=0
!p.linestyle=0

 cgcolorbar,range=range,annotatecolor='white',position=[0.03,0.07,0.57,0.085],charpercent=0.6,title=cbtitle,tickinterval=sigfig(tickinterval,2)


cgplots,[xsize-0.05*xsize],[0.05*ysize],color='white
cgplots,[xsize-0.05*xsize],[0.17*ysize],/continue,color='white'


cgplots,[xsize-0.05*xsize],[0.05*ysize],color='white
cgplots,[xsize-0.17*xsize],[0.05*ysize],/continue,color='white'

cgtext,xsize-0.2*xsize,ysize*0.04,'E',charthick=3,charsize=1.2,color=plotcolor

cgtext,xsize-0.06*xsize,ysize*0.18,'N',charthick=3,charsize=1.2,color=plotcolor

if keyword_set(platescale) then cgtext,xsize-0.13*xsize,ysize*0.02,string(sigfig(0.12*xsize*platescale,2))+'"',charthick=3,charsize=1.2,color=plotcolor

;end annotations

cgps_close,/png

cgcleanup

if do_2d then begin

cgps_open,location+name+ '_SNR_map.ps'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,xsize,ysize,/device

 ;  image = cgdemodata(18)
   cgImage,result,stretch=stretch,clip=clip_level,xrange=[0,xsize],yrange=[0,ysize],/keep_aspect_ratio,minvalue=-5.,maxvalue=10.,mean=meanvalue;,position=[0,0,1,1]



;begin annotations




!p.psym=0
!p.linestyle=0
  cgcolorbar,range=[-5,10],annotatecolor='white',position=[0.03,0.07,0.57,0.085],charpercent=0.6,title='SNR',tickinterval=3.


cgplots,[xsize-0.05*xsize],[0.05*ysize],color='white
cgplots,[xsize-0.05*xsize],[0.17*ysize],/continue,color='white'


cgplots,[xsize-0.05*xsize],[0.05*ysize],color='white
cgplots,[xsize-0.17*xsize],[0.05*ysize],/continue,color='white'

cgtext,xsize-0.2*xsize,ysize*0.04,'E',charthick=3,charsize=1.2,color=plotcolor

cgtext,xsize-0.06*xsize,ysize*0.18,'N',charthick=3,charsize=1.2,color=plotcolor

if keyword_set(platescale) then cgtext,xsize-0.13*xsize,ysize*0.02,string(sigfig(0.12*xsize*platescale,2))+'"',charthick=3,charsize=1.2,color=plotcolor

;end annotations

cgps_close,/png



endif ;2d if

cgps_open,location+name+ string(filter) +  '_source_map.ps'


   cgLoadCT, 1;, NColors=100,BOTTOM=100;,BOTTOM=100
	cgDisplay,xsize,ysize,/device

 ;  image = cgdemodata(18)
  ; cgImage,disp_frame,stretch=stretch,clip=clip_level,xrange=[0,xsize],yrange=[0,ysize],/keep_aspect_ratio,minvalue=minvalue,maxvalue=maxvalue,mean=meanvalue;,position=[0,0,1,1]


  if do_2d then cgImage,result,stretch=stretch,clip=clip_level,xrange=[0,xsize],yrange=[0,ysize],/keep_aspect_ratio,minvalue=-5,maxvalue=10,mean=meanvalue;,position=[0,0,1,1] else

if not do_2d then   cgImage,disp_frame,stretch=stretch,clip=clip_level,xrange=[0,xsize],yrange=[0,ysize],/keep_aspect_ratio,minvalue=minvalue,maxvalue=maxvalue,mean=meanvalue;,position=[0,0,1,1]


;begin annotations

if do_2d then cgcolorbar,range=[-5,10],annotatecolor='white',position=[0.03,0.07,0.57,0.085],charpercent=0.6,title='SNR',tickinterval=3.


if do_phot then range=[0,maxcon] else range=[minvalue,maxvalue]
if do_phot then cbtitle='Contrast' else cbtitle='ADU/Beam'
if do_phot then tickinterval=maxcon/5. else tickinterval=(maxvalue-minvalue)/5.
 


!p.psym=0
!p.linestyle=0


if not do_2d then  cgcolorbar,range=range,annotatecolor='white',position=[0.03,0.07,0.57,0.085],charpercent=0.6,title=cbtitle,tickinterval=sigfig(tickinterval,2)


cgplots,[xsize-0.05*xsize],[0.05*ysize],color='white
cgplots,[xsize-0.05*xsize],[0.17*ysize],/continue,color='white'


cgplots,[xsize-0.05*xsize],[0.05*ysize],color='white
cgplots,[xsize-0.17*xsize],[0.05*ysize],/continue,color='white'

cgtext,xsize-0.2*xsize,ysize*0.04,'E',charthick=3,charsize=1.2,color=plotcolor

cgtext,xsize-0.06*xsize,ysize*0.18,'N',charthick=3,charsize=1.2,color=plotcolor

if keyword_set(platescale) then cgtext,xsize-0.13*xsize,ysize*0.02,string(sigfig(0.12*xsize*platescale,2))+'"',charthick=3,charsize=1.2,color=plotcolor

;end annotations

;if where(sources_x eq 0.) eq -1 then begin

if n_elements(sources_x) eq 1 then if sources_x eq -1 then sources_x=[]
if n_elements(sources_y) eq 1 then if sources_y eq -1 then sources_y=[]
if n_elements(sources_contrast) eq 1 then if sources_contrast eq -1 then sources_contrast=[]


for jj=0, n_elements(sources_x)-1 do begin

			cgloadct
			PlotSym, 0, 2./(float(xsize)/1024.),thick=3
			!p.psym=8

			if sources_x[jj] ne 0. and sources_y[jj] ne 0. then begin

			cgplots,sources_x[jj]+0.5,sources_y[jj]+0.5,COLOR=plotcolor;,/device
			cgplots,sources_x[jj]+0.5,sources_y[jj]+0.5,/continue,COLOR=plotcolor;,/device


			cgtext,sources_x[jj]-35,sources_y[jj]+15,String(jj+1),charthick=3,charsize=1.2,color=plotcolor
			endif

endfor




;print, sources_contrast

;if sources_contrast ne !NULL then begin

if sources_contrast ne !NULL then if do_phot then sources_contrast=sources_contrast*correction_factor
if do_phot then tests=tests*correction_factor

;astrometry

if sources_contrast ne !NULL then if keyword_set(platescale) then sources_rho=sqrt( (sources_x-(xsize/2.))^2. + (sources_y-(ysize/2.))^2. ) * platescale else sources_rho=sqrt( (sources_x-(xsize/2.))^2. + (sources_y-(ysize/2.))^2. )
;print,  atan( (sources_y-(ysize/2.))/(sources_x-(xsize/2.)) )/!DTOR;+270.
;hak
if sources_contrast ne !NULL then sources_theta=atan( (sources_y-(ysize/2.)),(sources_x-(xsize/2.)) )/!DTOR+270.; mod 360.
;if sources_contrast ne !NULL then  sources_theta[where((sources_x-(xsize/2.)) gt 0)]=sources_theta[where((sources_x-(xsize/2.)) gt 0)]+180.
;if sources_contrast ne !NULL then sources_theta=atan( (sources_y-(ysize/2.)),(sources_x-(xsize/2.)) )/!DTOR+270.; mod 360.

;if sources_contrast ne !NULL then sources_theta[where(sources_x-(xsize/2.) lt 0)]=sources_theta[where(sources_x-(xsize/2.) lt 0)]+180.

if sources_contrast ne !NULL then sources_theta=sources_theta mod 360.

printf,1, '# ------[ Source list for '+name+']------'
printf,1, '# --------- Minimum significance = ',sigma
printf,1, '#-----------Corr_thresh = ',ct 
if keyword_set(platescale) then printf,1, '# --------- Plate scale ("/px) = ',platescale
;printf,1, '# --------- Wavelength (m) = ',wavelength
;printf,1, '# --------- Diameter (m) = ',diameter
printf,1, '# --------- FWHM (px) = ',FWHM
if do_phot then printf,1, '# --------- Throughput correction factor = ',correction_factor

printf,1, '# --------------------------------------------------------------'
if keyword_set(platescale) then begin 
	if do_phot then printf,1, '# Source 	x		y	sigma     contrast	  rho (")	theta' else printf,1, 'Source #	x		y	sigma	     rho (")	   theta'
endif
if not keyword_set(platescale) then begin
 if do_phot then printf,1, '# Source 	x		y	sigma     contrast	  rho (px)	theta' else printf,1, 'Source #	x		y	sigma	     rho (px)	   theta'

endif
printf,1, '# --------------------------------------------------------------'

if sources_contrast ne !NULL then for nn=0,n_elements(sources_x)-1 do $
	if do_phot then printf,1,nn+1,sources_x[nn],sources_y[nn],sources_sigma[nn], sources_contrast[nn],sources_rho[nn],sources_theta[nn]$
		 else printf,1,nn+1,sources_x[nn],sources_y[nn],sources_sigma[nn],sources_rho[nn],sources_theta[nn]

printf,1, '# --------------------------------------------------------------'
if sources_contrast ne !NULL then if do_phot then printf,1, '# --- Maximum source contrast',max(sources_contrast)
if sources_contrast ne !NULL then if do_phot then printf,1, '# --- Minimum source contrast',min(sources_contrast)
if sources_contrast ne !NULL then if do_phot then printf,1, '# --- Median source contrast',median(sources_contrast)
if sources_contrast ne !NULL then if do_phot then printf,1, '# --- Mean source contrast',mean(sources_contrast)
if n_elements(sources_contrast) gt 20 and do_phot then mmm,sources_contrast,mode
if sources_contrast ne !NULL then if n_elements(sources_contrast) gt 20 and do_phot then printf,1, '# --- Mode source contrast',mode
if sources_contrast ne !NULL then if do_phot then printf,1, '# --- Stddev source contrast',stddev(sources_contrast)
if sources_contrast ne !NULL then if do_phot then printf,1, '# --------------------------------------------------------------'
if sources_contrast ne !NULL then printf,1, '# --- Maximum source significance: ',max(sources_sigma)
if sources_contrast ne !NULL then printf,1, '# --- Minimum source significance: ',min(sources_sigma)
if sources_contrast ne !NULL then printf,1, '# --- Mean source significance: ',mean(sources_sigma)

save,filename=location+name+ '_sources.sav',sources_x,sources_y,sources_sigma,sources_contrast,sources_rho,sources_theta


;hak
	;cgcleanup

;tvlct,255,255,255,0
;!p.color=0
cgps_close,/png
cgps_open,location+name+ '_radial_sensitivity.ps',/landscape
	;window,1,xs=xsize,ys=ysize
!p.psym=0

if keyword_set(platescale) then xtitle='Radius (Arcsec)' else xtitle= 'Radius (px)'

sensitivity=[]
;for ii=0,n_elements(tests)-1 do sensitivity=[sensitivity,5.*(tests[ii]-testmeans[ii])/ref_flux]
for ii=0,n_elements(tests)-1 do sensitivity=[sensitivity,(sigma-1.)*(tests[ii])/double(ref_flux)] ;sigma - 1 to fix the bug of +5sigma points appearing below the line. The cause of the bug has not been found as of 20180101
if keyword_set(platescale) then xrange=[0,(inrad*platescale)+(n_elements(sensitivity)*platescale)] else xrange=[0,n_elements(sensitivity)]


if do_phot then yrange=[1E-7,1E-0] else yrange=[min(sensitivity,/nan)/2.,1.5*max(sensitivity,/nan)]
if fwhm gt 5 then yrange=[1E-5,1E-2]
if do_phot then ytitle=textoidl(string(sigfig(5.0,2))+'-\sigma Contrast') else ytitle = textoidl(string(sigfig(sigma,2))+'-\sigma Sensitivity (ADU/Beam)')

if not keyword_set(platescale) then platescale=1.


plot, (inrad*platescale)+(findgen(n_elements(sensitivity))*platescale),sensitivity,charsize=2,ytitle=ytitle,xtitle=xtitle,/ylog,xrange=xrange,yrange=yrange
!p.psym=4
if sources_contrast ne !NULL then oplot, sqrt( (sources_x-(xsize/2.))^2. + (sources_y-(ysize/2.))^2. )*platescale, sources_contrast


if sources ne !NULL then writefits,location+name+ '_sources.fits',sources

find_correction_value=0	;useful for determining the radially dependent correction factor
true_contrast=(3.0E-5)	;known contrast of the majority of sources if the field of view
			;if there are a lot of background stars, inject a lot of sources

;if find_correction_value then begin
;	yerr=fltarr(n_elements(sources_rho))
;	yerr=(5./sources_sigma)*stddev(sources_contrast)
;	sources_rho_px=sources_rho/platescale
 ;  coeffs = POWERLAW_FIT(sources_rho_px,sources_contrast,errors=yerr)
	

;print, coeffs
;!p.psym=0
;oplot,findgen(xsize/2.),(findgen(xsize/2.)/coeffs[0])^(coeffs[1])

;correction=true_contrast/((findgen(xsize/2.)/coeffs[0])^(coeffs[1]))


;writefits,location+name+'_radial_throughput_correction.fits',correction
;endif

cgps_close,/png


if platescale lt 1 then begin
cgps_open,location+name+ '_radial_sensitivity_2.0_arcsec.ps',/landscape
	;window,1,xs=xsize,ys=ysize
!p.psym=0

if keyword_set(platescale) then xtitle='Radius (Arcsec)' else xtitle= 'Radius (px)'

;sensitivity=[]
;for ii=0,n_elements(tests)-1 do sensitivity=[sensitivity,sigma*tests[ii]/ref_flux]
;if keyword_set(platescale) then xrange=[0,n_elements(sensitivity)*platescale] else xrange=[0,n_elements(sensitivity)]
xrange=[0,2]

if do_phot then yrange=[1E-7,1E-0] else yrange=[min(sensitivity,/nan)/2.,1.5*max(sensitivity,/nan)]
if do_phot then ytitle=textoidl(string(sigfig(5.0,2))+'-\sigma Contrast') else ytitle = textoidl(string(sigfig(sigma,2))+'-\sigma Sensitivity (ADU/Beam)')

if not keyword_set(platescale) then platescale=1.


plot, (inrad*platescale)+(findgen(n_elements(sensitivity))*platescale),sensitivity,charsize=2,ytitle=ytitle,xtitle=xtitle,/ylog,xrange=xrange,yrange=yrange
!p.psym=4
if sources_contrast ne !NULL then oplot, sqrt( (sources_x-(xsize/2.))^2. + (sources_y-(ysize/2.))^2. )*platescale, sources_contrast


cgps_close,/png;,width=1500
endif

;endif ;sources ne null if

;endif ;sources_x


;	cgcleanup

;close, 1; Commented out because of writing to txt error





end
