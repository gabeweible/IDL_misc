pro kklip_loop_snr, coadd, kklip_arr=kklip_arr, angsep_arr=angsep_arr, bin_arr=bin_arr, extra=extra
;+
; NAME:
;       KKLIP_LOOP_SNR

; PURPOSE:
;       Calculate SNR for reductions performed over various parameters (e.g., K_KLIP)

;       WRITTEN, 2025 Gabriel Weible
;-

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)

output_path = '/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/Alcor/macbook_' +$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)

restore, output_path+'/kklip_loop_arrays.sav'
fwhm=13.0773; M-band]
	
bin_type_arr = ['median', 'mean']
combine_type_arr = ['median']
;filter_arr = 100*randomu(seed, 5);[30., 40., 50., 60., 70., 80., 90.]
n_ang_arr = [1,2,4]
;if not keyword_set(bin_arr) then szbin=1 & bin=szbin
bin_arr = [4];I should probably try lower coadds as well!
pa_off=0.
;lambda_over_D = 11.07 ; px
coadd=50
wr = 34-7


; NOTE: LAMBDA/D ~ 11.07 PX, SHOULD HAVE 2PI*(R/11.07PX) RESOLUTION ELEMENTS: 1 FOR ALCOR B AND THE
; OTHERS FOR NOISE COMPUTATION
; ALCOR B AT ~10.45 PX RADIUS SO 0.945 LAMBDA/D SO 6 RESOLUTION ELEMENTS, 1 FOR ALCOR B, 5 FOR NOISE

results = []
files = []

pre_inj=0
neg_inj=0
uncert=0
klip=1
use_gauss=0
fs=0
nod='dx_only'
cds='endcds'
pre_clean=0
adi=0
do_annmode=1

anglemax = 360.
nrings=1
pxscale_dx = 0.010610; 
pxscale = pxscale_dx
xp = 176.6-1.; back to no pa_off;189.163; approx. x-pixel position of Alcor B (w/pa_off=22.)
yp = 207.3-1.; back to no pa_off;198.093 ; approx. y-pixel position of Alcor B (w/pa_off=22.)

half_cropped_sz = 190.
full_frame = 0

if full_frame eq 1 then wr=fix(half_cropped_sz/float(nrings)); width of ring for KLIP

corr_thresh = 0.9996
obj = 'Alcor'
cube_folder = output_path

best_sig = 0.
best_sig2 = 0.
best_sig3 = 0.
best_sig_file2 = ''
best_sig_file3 = ''
best_sig_file = ''
best_signal = 0.
best_signal2 = 0.
best_signal3= 0.
best_signal_file = ''
best_signal_file2 = ''
best_signal_file3 = ''

; Loop through all pixels in the image
foreach szbin, bin_arr do begin
		foreach n_ang, n_ang_arr do begin
		foreach angsep, angsep_arr do begin
		foreach filter, filter_arr do begin
		foreach bin_type, bin_type_arr do begin
		foreach combine_type, combine_type_arr do begin
		foreach k_klip, kklip_arr do begin
		
			print, 'running Alcor SNR calc. for k_klip = ', k_klip
			
			suffix=strcompress(reform('_'+string(k_klip)+'k_'+string(sigfig(angsep,4))+$
				'as_'+string(sigfig(anglemax,2))+'am_'+string(sigfig(nrings,2))+'rings_'+$
				string(wr)+'wr_'+String(n_ang)+'nang_'+string(sigfig(filter,3))+'filter_'+$
				string(sigfig(szbin,2))+'bin'+string(corr_thresh)+'corrthresh'), /remove_all)
	 			
	 		;read image
	 		
	 		;Alcor_left_klip_8k_0.2218as_360am_1.0rings_190wr_1nang_13.filter_6.0bin0.997000corrthresh_bin_6._type_mean_comb_type_nwadi_k_klip_8.00_angsep_0.2218_angmax_360._nrings_1.0_nang_1.0_neg_inj_0
	 		
	 		file = strcompress(cube_folder + '/combined/' + obj + '_right_klip' + suffix + '_bin_' +$
   				string(sigfig(szbin,1)) + '_type_' + bin_type + '_comb_type_'+combine_type+'_k_klip_'+$
   				string(sigfig(k_klip,3))+'_angsep_'+string(sigfig(angsep,4))+'_angmax_'+$
   				string(sigfig(anglemax,3))+'_nrings_'+string(sigfig(nrings, 2))+'_nang_'+$
   				string(sigfig(n_ang,2)) + '_neg_inj_' + string(neg_inj) +   '.fits', /rem)
   				
			image=readfits(file)
		
		; remove zeros
		image[where(image eq 0)] = !values.f_nan

		;what is rho, theta on the image?
		xs = (size(image))[1] & ys = (size(image))[2]
		xhs = xs/2. & yhs = ys/2. ;half-sizes
		print, xhs, yhs
		
		print, 'Xp,Yp = ',xp,yp
		
		xc=xhs & yc=yhs  & xsize=xs & ysize=ys
		print, xhs, yhs
		
		; NOISE IMAGE SETUP
		nimage = image ;make a copy
		boxehs = fwhm/2. ;source exclusion half-size
		; Exclude injected source for the 'noise' image
		nimage[xp-boxehs:xp+boxehs-1, yp-boxehs:yp+boxehs-1] = !values.f_nan
		
		; get rid of really high or low values (outside 10 sigma from the mean)
		nantr=10
		nanmask=where(nimage lt (mean(nimage,/nan,/double)-nantr*stddev(nimage,/nan,/double) ) or nimage gt (mean(nimage,/nan,/double)+nantr*stddev(nimage,/nan,/double) ))
		image[nanmask]=!values.f_nan
		nimage[nanmask]=!values.f_nan
		
		print, fwhm; & hak

		;measure signal in an aperture
		aper_rad=fwhm/2.;1.      ;fwhm ;& skyradi=4. & skyrado=4.
		zimage=image
		zimage[where(finite(zimage) ne 1)]=0.
		
		; fwhm-diameter aperture for signal, but we'll do 1/2 fwhm-diamter for "significance"
		aper, zimage, xp,yp, signal,err,sky,skyerr,1,aper_rad,[0,0],[-99E99,99E99],/flux,/exact,/silent,SETSKYVAL = 0,/nan;signal
		
		;copied from find_sources
		;we'll calculate flux in each aperture (out of ntheta total)
		print, 'xp:', xp
		print, 'yp', yp
		rho_px = ((xhs - xp)^2. + (yhs - yp)^2.)^(1./2.); radius from center of image in pixels
		theta = acos((xp-xhs) / rho_px)
		print, 'theta of target:', theta
		;ntheta = 5;4 for noise (fwhm - +1 if lambda/D)
		ntheta = float(fix(2.*!PI*float(rho_px)/FWHM))
		;ntheta=float( round( ( 2.*!PI*float(rho_px)/(2.0*aper_rad) ) )) ;about one angular bin per resolution element
		print, 'ntheta:', ntheta
		nzimage=nimage
		nzimage[where(finite(nzimage) ne 1)]=0.; get rid of nans

		if ntheta lt 6 then begin stitheta=1 & stendtheta=ntheta-1 & endif else begin & stitheta=2 & stendtheta=ntheta-2 & endelse
		;stitheta = 1 & stendtheta = ntheta-1
		for itheta=stitheta,stendtheta do begin
			ttheta=theta+(float(itheta)/float(ntheta))*2.*!PI;+(!PI/31.)
			xtest=xhs+(float(rho_px)*cos(ttheta)) & ytest=yhs+(float(rho_px)*sin(ttheta))

			aper, nzimage, xtest, ytest, flux,fluxerr,sky,skyerr,1,aper_rad,[0,0],[-99E99,99E99],/flux,/exact,/silent,SETSKYVAL = 0,/nan; noise at itheta aperture

			; append to arrays for all the tests
			if itheta eq stitheta then fluxarr=flux else fluxarr=[fluxarr,flux]
			if itheta eq stitheta then xarr=xtest else xarr=[xarr,xtest]
			if itheta eq stitheta then yarr=ytest else yarr=[yarr,ytest]
		endfor ;;itheta for
		
		fluxarr[(where(fluxarr eq 0.))]=!values.f_nan
		signal=signal-mean(fluxarr,/nan,/double); signal is a difference from mean flux at this radius
		noise=stddev(fluxarr,/nan,/double)*sqrt(1.+(1./n_elements(fluxarr[where(finite(fluxarr) ne 0)]))) ;small sample statistics correction: Mawet+2014)
		print, 'signal', signal
		result=signal/noise; "Real" SNR from Mawet+ 2014, maybe should be called a "significance" (corresponds to Gaussian sigma statistics)
		print, 'Source recovered with significance = ',result

		; keep track of the best signal on Alcor B (running)
		if signal gt best_signal3 then begin
			if signal gt best_signal2 then begin
				if signal gt best_signal then begin; better than 1st
					best_signal3 = best_signal2
					best_signal2 = best_signal
					best_signal = signal
			
					best_signal_file3 = best_signal_file2
					best_signal_file2 = best_signal_file
					best_signal_file = file
				endif else begin; better than 2nd, worse than 1st
					best_signal3 = best_signal2
					best_signal2 = signal
			
					best_signal_file3 = best_signal_file2
					best_signal_file2 = file
				endelse; better than 2nd else
			endif else begin; better than 3rd, but not 2nd
				best_signal3 = signal
				best_signal_file3 = file
			endelse; better than 3rd
		endif; better than top3 end
				
		; keep track of the best significance at Alcor B (running)
		results = [results,result] ;store results
		if result gt best_sig3 then begin
			if result gt best_sig2 then begin
				if result gt best_sig then begin; better than 1st
					best_sig3 = best_sig2
					best_sig2 = best_sig
					best_sig = result
			
					best_sig_file3 = best_sig_file2
					best_sig_file2 = best_sig_file
					best_sig_file = file
				endif else begin; better than 2nd, worse than 1st
					best_sig3 = best_sig2
					best_sig2 = result
			
					best_sig_file3 = best_sig_file2
					best_sig_file2 = file
				endelse; better than 2nd else
			endif else begin; better than 3rd, but not 2nd
				best_sig3 = result
				best_sig_file3 = file
			endelse; better than 3rd
		endif; better than top3 end
		
		print, newline, 'best significance:', best_sig
		print, 'best sig file:', best_sig_file
		print, 'best signal:', best_signal
		print, 'best signal file:', best_signal_file, newline
		
		print,'best significance2:', best_sig2
		print, 'best sig file2:', best_sig_file2
		print, 'best signal2:', best_signal2
		print, 'best signal file2:', best_signal_file2, newline
		
		print,'best significance3:', best_sig3
		print, 'best sig file3:', best_sig_file3
		print, 'best signal3:', best_signal3
		print, 'best signal file3:', best_signal_file3, newline
		
		files = [files, file]
		
		save,filename=strcompress(output_path + '/Alcor_contrast_results.sav',/rem),$
			results, files
			
		; Write curves to a CSV.
		WRITE_CSV, strcompress(output_path + '/Alcor_contrast_output.csv',/rem), results
		
	 	
	 	endforeach; k_klip loop
	 	endforeach; combine_type
	 	endforeach; bin_type loop
		endforeach; filter loop
	 	endforeach; angsep loop
	 	endforeach; nangarr
endforeach;szbin loop

end