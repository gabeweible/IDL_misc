pro kklip_loop8, kklip_arr=kklip_arr, angsep_arr=angsep_arr, bin_arr=bin_arr
;alcor_pipeline, pre_inj=1,pre_clean=1,coadd=25,klip=1,adi=0,fs=0,nod='dx_only',cds='endcds',k_klip=5,do_annmode=1,szbin=4,angsep=1,bin_type='median',filter=25.,n_ang=1,pa_off=0,combine_type='median'

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

fwhm=13.0773; M-band
;lambda_over_D = 11.07 ; px
coadd=25
pa_off=0.

; this should roughly go fastest to slowest, except for the
; loop over kklip (which I think is important)
; I probably should loop over angsep deeper into the nested for loops
; so that this relatively important parameter may be for quickly varied, as well

kklip_arr=[26,27,28,29,30,31];[4,5,6,7];[4,5,6,7];new array from initial results;[8,9]; partially complete;[2, 3, 4, 5, 6, 7, 8, 9]
angsep_arr=[0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51];[0.05, 0.06, 0.07, 0.7, 0.8, 0.9];[0.25/fwhm, 0.5/fwhm, 0.75/fwhm, 1.0/fwhm, 2.0/fwhm, 3.0/fwhm, 3.25/fwhm, 3.5/fwhm, 3.75/fwhm, 4.0/fwhm, 5.0/fwhm]

bin_type_arr = ['median']
combine_type_arr = ['median']
filter_arr = [29,33,37,41,45,49,53,57,61,65,69,73,77,81];[30., 40., 50., 60., 70., 80., 90.]
n_ang_arr = [1]
bin_arr = [4];I should probably try lower coadds as well!

output_path = '/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/Alcor/macbook_' +$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)

save, filename=output_path+'/kklip_loop_arrays8.sav',kklip_arr, angsep_arr, filter_arr,$
	bin_type_arr, combine_type_arr, n_ang_arr, bin_arr

; Loop through all pixels in the image
;foreach sky_type, sky_type_arr do begin
foreach szbin, bin_arr do begin
		foreach bin_type, bin_type_arr do begin
		foreach combine_type, combine_type_arr do begin
		foreach n_ang, n_ang_arr do begin
		foreach filter, filter_arr do begin
		foreach angsep, angsep_arr do begin
		foreach k_klip, kklip_arr do begin
		
			print, newline, 'running Alcor pipeline in loop for:',newline
			;print, newline, 'sky_type = ', sky_type
			print, 'szbin = ', string(szbin)
			print, 'n_ang = ', string(n_ang)
			print, 'angsep = ', string(angsep)
			print, 'filter= ', string(filter)
			print, 'bin_type: ', bin_type
			print, 'combine_type:', combine_type
			print, 'k_klip = ', string(k_klip), newline
			
				alcor_pipeline, pre_inj=0, neg_inj=0, coadd=coadd, uncert=0, klip=1,$
	 			use_gauss=0, fs=0, nod='dx_only', cds='endcds', pre_clean=0,$
	 			k_klip=k_klip, adi=0, do_annmode=1, szbin=szbin, angsep=angsep,$
	 			bin_type=bin_type, filter=filter, n_ang=n_ang, pa_off=pa_off,$
	 			combine_type=combine_type
	 			
	 	endforeach; k_klip loop
	 	endforeach; angsep loop
		endforeach; filter loop
	 	endforeach; nangarr
	 	endforeach; combine_type
	 	endforeach; bin_type loop
endforeach;szbin loop
;endforeach; sky_type loop

end