pro kklip_loop2, kklip_arr=kklip_arr, angsep_arr=angsep_arr, bin_arr=bin_arr
;+
; NAME:
;       MASK_RESIDUALS

; PURPOSE:
;       Reduce PSF-subtracted residuals in an the center of an image to zero,
;       within a radius

; CALLING SEQUENCE:
;       MASK_RESIDUALS, image, radius, output_file

; INPUTS:
;       image = a single 2D FITS image
; 
;       radius = radius in pixels to mask in a circle (15 good for HII 1348B)

;       output_file = string output FITS to write the masked image to
;       
; RESTRICTIONS:
;       (1) image must be a string path to the FITS image, radius may be a float
;           or an integer.
;
; EXAMPLE:
;       Mask the center of an image within a radius of 10 pixels:
;       MASK_RESIDUALS, '/path/to/image.fits', 10, '/path/to/masked/image.fits'
;
; EXTERNAL PROCEDURES USED:
;       READFITS, WRITEFITS
;
; MODIFICATION HISTORY:
;       WRITTEN, 2024 Gabriel Weible
;-
; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)

kklip_arr=[71,72,73,74,76,77,78,79,81,82];[4,5,6,7];[4,5,6,7];new array from initial results;[8,9]; partially complete;[2, 3, 4, 5, 6, 7, 8, 9]
fwhm=13.0773; M-band
angsep_arr=0.2+0.6*randomu(seed, 10);[0.05, 0.06, 0.07, 0.7, 0.8, 0.9];[0.25/fwhm, 0.5/fwhm, 0.75/fwhm, 1.0/fwhm, 2.0/fwhm, 3.0/fwhm, 3.25/fwhm, 3.5/fwhm, 3.75/fwhm, 4.0/fwhm, 5.0/fwhm]

	
bin_type_arr = ['mean']
combine_type_arr = ['median', 'nwadi']
filter_arr = 30 + 80*randomu(seed, 10);[30., 40., 50., 60., 70., 80., 90.]
n_ang_arr = [2,4]
;if not keyword_set(bin_arr) then szbin=1 & bin=szbin
bin_arr = [2,1];I should probably try lower coadds as well!
pa_off=0.
;lambda_over_D = 11.07 ; px
coadd=50

output_path = '/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/Alcor/macbook_' +$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)

save, filename=output_path+'/kklip_loop_arrays2.sav',kklip_arr, angsep_arr, filter_arr

; Loop through all pixels in the image
foreach szbin, bin_arr do begin
		foreach n_ang, n_ang_arr do begin
		foreach angsep, angsep_arr do begin
		foreach filter, filter_arr do begin
		foreach bin_type, bin_type_arr do begin
		foreach combine_type, combine_type_arr do begin
		foreach k_klip, kklip_arr do begin


	;print, 'running HII 1348B pipeline for k_klip = ', k_klip
	;hii1348_pipeline, pre_inj=0, neg_inj=0, coadd=25, uncert=0, klip=1,$
		;use_gauss=0, extra='kklip_'+strcompress(k_klip, /r), k_klip=k_klip,$
		;nod='total'
		
			print, 'running Alcor pipeline for k_klip = ', k_klip
				alcor_pipeline, pre_inj=0, neg_inj=0, coadd=coadd, uncert=0, klip=1,$
	 			use_gauss=0, fs=0, nod='dx_only', cds='endcds', pre_clean=0,$
	 			k_klip=k_klip, adi=0, do_annmode=1, szbin=szbin, angsep=angsep, bin_type=bin_type,$
	 			filter=filter, n_ang=n_ang, pa_off=pa_off, combine_type=combine_type
	 	endforeach; k_klip loop
	 	endforeach; combine_type
	 	endforeach; bin_type loop
		endforeach; filter loop
	 	endforeach; angsep loop
	 	endforeach; nangarr
endforeach;szbin loop

end