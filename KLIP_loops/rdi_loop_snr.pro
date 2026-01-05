pro rdi_loop_snr, coadd, kklip_arr=kklip_arr, angsep_arr=angsep_arr,$
	bin_arr=bin_arr, extra=extra, filename=filename
;+
; NAME:
;       KKLIP_LOOP_SNR

; PURPOSE:
;       Calculate SNR for reductions performed over various parameters (e.g., K_KLIP)

;       WRITTEN, 2025 Gabriel Weible
;-

compile_opt IDL2
newline = string(10B)

output_path = '/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/Alcor/macbook_' +$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)

if not keyword_set(filename) then begin
	; 2-digit k_klip
	files = FILE_SEARCH(output_path, '/Alcor_k_klip_[0-9][0-9].0_comb_type_*_combined.fits', COUNT=filecount)
	print, 'Found ', filecount, 'right KLIP FITS files!'
	
	; one-digit k_klip
	files2 = FILE_SEARCH(output_path, '/Alcor_k_klip_[0-9].00_comb_type_*_combined.fits', COUNT=filecount2)
	print, 'Found ', filecount2, 'right KLIP FITS files!'
	
	; 3-digit k_klip
	;files3 = FILE_SEARCH(output_path, '/Alcor_k_klip_[0-9][0-9][0-9]._comb_type_*_combined.fits', COUNT=filecount3)
	;print, 'Found ', filecount3, 'right KLIP FITS files!'
	
	; 4-digit k_klip
	;files4 = FILE_SEARCH(output_path, '/Alcor_right_rdi_[0-9][0-9][0-9][0-9]k*_neg_inj_0.fits', COUNT=filecount4)
	;print, 'Found ', filecount4, 'right KLIP FITS files!'
	
	all_files = [files, files2];, files3];, files4]
endif else begin
	all_files = [filename]
endelse

avg_results = []

fwhm=13.0773; M-band]
xp_c = 179.25-1.; back to no pa_off;189.163; approx. x-pixel position of Alcor B (w/pa_off=22.)
yp_c = 208.25-1.; back to no pa_off;198.093 ; approx. y-pixel position of Alcor B (w/pa_off=22.)
xp_arr = [xp_c-1., xp_c-0.5, xp_c, xp_c+0.5, xp_c+1.]
yp_arr = [yp_c-1., yp_c-0.5, yp_c, yp_c+0.5, yp_c+1.]

; Note: radius is ~ 21.9 px for Alcor B = 1.98 lambda/D = 1.675 fwhm (!!!)
; Note: r = 232.465 mas
; Alcor is at 25.056 pc from Hipparcos, so this is:
; 5.825 au projected separation!!!! Jupiter is at 5.20 au !!!!
half_cropped_sz = 190.
full_frame = 0

if full_frame eq 1 then wr=fix(half_cropped_sz/float(nrings)); width of ring for KLIP

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

foreach file,all_files do begin
	;print, newline, 'reading: ', file, newline
   				
	image=readfits(file)
	
	; remove zeros
	image[where(image eq 0)] = !values.f_nan

	;what is rho, theta on the image?
	xs = (size(image))[1] & ys = (size(image))[2]
	xhs = xs/2. & yhs = ys/2. ;half-sizes
	
	
	test_grid_signals = []
	test_grid_results = []
	foreach yp, yp_arr do begin
	foreach xp, xp_arr do begin
		
		xc=xhs & yc=yhs  & xsize=xs & ysize=ys
		;print, xhs, yhs
		
		; NOISE IMAGE SETUP
		nimage = image ;make a copy
		boxehs = fwhm/2. ;source exclusion half-size
		; Exclude the expected source for the 'noise' image
		nimage[xp-boxehs:xp+boxehs-1, yp-boxehs:yp+boxehs-1] = !values.f_nan
		
		; get rid of really high or low values (outside 10 sigma from the mean)
		nantr=10
		nanmask=where(nimage lt (mean(nimage,/nan,/double)-nantr*stddev(nimage,/nan,/double) ) or nimage gt (mean(nimage,/nan,/double)+nantr*stddev(nimage,/nan,/double) ))
		image[nanmask]=!values.f_nan
		nimage[nanmask]=!values.f_nan
	
		;measure signal in an aperture
		aper_rad=fwhm/2.;1.      ;fwhm ;& skyradi=4. & skyrado=4.
		zimage=image
		zimage[where(finite(zimage) ne 1)]=0.
		
		; fwhm-diameter aperture for signal, but we'll do 1/2 fwhm-diamter for "significance"
		aper, zimage, xp,yp, signal,err,sky,skyerr,1,aper_rad,[0,0],[-99E99,99E99],/flux,/exact,/silent,SETSKYVAL = 0,/nan;signal
		
		;copied from find_sources
		;we'll calculate flux in each aperture (out of ntheta total)
		;print, 'xp:', xp
		;print, 'yp', yp
		rho_px = ((xhs - xp)^2. + (yhs - yp)^2.)^(1./2.); radius from center of image in pixels
		theta = acos((xp-xhs) / rho_px)
		;print, 'theta of target:', theta
		;ntheta = 5;4 for noise (fwhm - +1 if lambda/D)
		ntheta = float(fix(2.*!PI*float(rho_px)/FWHM))
		;ntheta=float( round( ( 2.*!PI*float(rho_px)/(2.0*aper_rad) ) )) ;about one angular bin per resolution element
		;print, 'ntheta:', ntheta
		nzimage=nimage
		nzimage[where(finite(nzimage) ne 1)]=0.; get rid of nans
	
		if ntheta lt 6 then begin
			stitheta=1
			stendtheta=ntheta-1
		endif else begin 
			stitheta=2
			stendtheta=ntheta-2
		endelse
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
		;print, 'signal', signal
		result=signal/noise; "Real" SNR from Mawet+ 2014, maybe should be called a "significance" (corresponds to Gaussian sigma statistics)
		print, 'Source recovered with significance = ', result
		
		test_grid_signals = [test_grid_signals, signal]
		test_grid_results = [test_grid_results, result]
	endforeach; loop over xp_arr
	endforeach; loop over yp_arr
	
	avg_signal = median(test_grid_signals, /double, /even);; mean(test_grid_signals, /double, /nan); max(test_grid_signals);
	avg_result = median(test_grid_results, /double, /even);; mean(test_grid_results, /double, /nan);  max(test_grid_signals);

	; keep track of the best signal on Alcor B (running)
	if avg_signal gt best_signal3 then begin
		if avg_signal gt best_signal2 then begin
			if avg_signal gt best_signal then begin; better than 1st
				best_signal3 = best_signal2
				best_signal2 = best_signal
				best_signal = avg_signal
		
				best_signal_file3 = best_signal_file2
				best_signal_file2 = best_signal_file
				best_signal_file = file
			endif else begin; better than 2nd, worse than 1st
				best_signal3 = best_signal2
				best_signal2 = avg_signal
		
				best_signal_file3 = best_signal_file2
				best_signal_file2 = file
			endelse; better than 2nd else
		endif else begin; better than 3rd, but not 2nd
			best_signal3 = avg_signal
			best_signal_file3 = file
		endelse; better than 3rd
	endif; better than top3 end
			
	; keep track of the best significance at Alcor B (running)
	avg_results = [avg_results, avg_result] ;store results
	
	if avg_result gt best_sig3 then begin
		if avg_result gt best_sig2 then begin
			if avg_result gt best_sig then begin; better than 1st
				best_sig3 = best_sig2
				best_sig2 = best_sig
				best_sig = avg_result
		
				best_sig_file3 = best_sig_file2
				best_sig_file2 = best_sig_file
				best_sig_file = file
			endif else begin; better than 2nd, worse than 1st
				best_sig3 = best_sig2
				best_sig2 = avg_result
		
				best_sig_file3 = best_sig_file2
				best_sig_file2 = file
			endelse; better than 2nd else
		endif else begin; better than 3rd, but not 2nd
			best_sig3 = avg_result
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
		
endforeach;files loop

; Write curves to a CSV.
WRITE_CSV, strcompress(output_path + '/Alcor_contrast_output_big_rdi.csv',/rem), avg_results, all_files
save,filename=strcompress(output_path + '/Alcor_contrast_results_big_rdi.sav',/rem),$
	avg_results, all_files

end