pro kklip_loop_adi_new2, kklip_arr=kklip_arr, angsep_arr=angsep_arr, bin_arr=bin_arr

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)

;extra = 'new'

fwhm=13.0773; M-band
;lambda_over_D = 11.07 ; px
coadd=20
pa_off=0.
;save_basis=1; start by saving (only necessary as filter, bin change)
save_basis=0; don't save the basis for "RDI" that's really klip...
do_dewarp=1; don't dewarp if testing RDI

bin_type_arr = ['mean']
combine_type_arr = ['median']

; hmmm...I might want to try pretty conservative values here to retain signal.
filter_arr = [0.0]
;lp_filter_arr = [0];1.0;fwhm*0.25; low-pass filter

bin_arr = [1];I should probably try lower coadds as well!

keep_number_arr = [1000, 1025, 1050, 1075, 1100, 1125, 1150, 1175, 1200]
	
;peak_thresh_arr = [0.9]
	
stddev_thresh_arr = [0.8, 0.85, 0.9, 0.95, 1]
	
normal=0

obj='Alcor'

output_path = '/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/Alcor/macbook_' +$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)

; start to loop through parameters, where saving a new basis and binned cubes for
; for each new sz bin and/or filter size

foreach stddev_thresh, stddev_thresh_arr do begin
pre_inj=1

;foreach peak_thresh, peak_thresh_arr do begin
;pre_inj=1

foreach keep_number, keep_number_arr do begin
pre_inj=1
	
foreach szbin, bin_arr do begin
	foreach bin_type, bin_type_arr do begin
		foreach combine_type, combine_type_arr do begin
				foreach filter, filter_arr do begin
					;foreach lp_filter, lp_filter_arr do begin

					print, newline, 'running Alcor pipeline in loop for:',newline
					print, 'szbin = ', string(szbin)
					print, 'filter= ', string(filter)
					print, 'bin_type: ', bin_type
					print, 'combine_type:', combine_type
					
					alcor_pipeline, pre_inj=pre_inj, neg_inj=0, coadd=coadd,$
						uncert=0, klip=0, use_gauss=0, fs=0, nod='dx_only',$
						pre_clean=0, adi=1, szbin=szbin, bin_type=bin_type,$
						filter=filter, combine_type=combine_type,do_rdi=0,$
						do_dewarp=do_dewarp,$
						keep_number=keep_number, sky_type='median',$
						do_create_basis=1, do_new_binning=1, $
						save_basis=0, pa_off=pa_off,$
						normal=normal, peak_clean=0, stddev_clean=1,$
						stddev_thresh=stddev_thresh
						
					pre_inj=0
					
					;endforeach; lp_filter
								
				endforeach; filter loop
		endforeach; combine_type
	endforeach; bin_type loop
endforeach;szbin loop

endforeach; corr_thresh loop
;endforeach; peak_thresh loop
endforeach; stddev_thresh loop

end