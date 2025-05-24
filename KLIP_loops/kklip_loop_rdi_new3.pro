pro kklip_loop_rdi_new3, kklip_arr=kklip_arr, angsep_arr=angsep_arr, bin_arr=bin_arr

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)

extra = 'new3'

lp_filter = 1; low-pass filter with 2-pixel width

fwhm=13.0773; M-band
;lambda_over_D = 11.07 ; px
coadd=20
pa_off=0.
;save_basis=1; start by saving (only necessary as filter, bin change)
save_basis=0; don't save the basis for "RDI" that's really klip...
do_dewarp=0; don't dewarp if testing RDI

; this should roughly go fastest to slowest, except for the
; loop over kklip (which I think is important)
; I probably should loop over angsep deeper into the nested for loops
; so that this relatively important parameter may be for quickly varied, as well

kklip_arr = [4,8,16,32,64,128]
angsep_arr=[0.0, 0.15, 0.3, 0.45, 0.6]

bin_type_arr = ['mean']
combine_type_arr = ['median']
filter_arr = [0.0];[30., 40., 50., 60., 70., 80., 90.]
n_ang_arr = [1]
bin_arr = [1];I should probably try lower coadds as well!

; this should be 2x FWHM in width after everything is actually combined
annmode_inout=[8, 34];
annmode_in = annmode_inout[0]
annmode_out = annmode_inout[1]
r_px = ((annmode_out + annmode_in) / 2.)

; 1-px of overlap. Then, I just take 2nd-from-inside-most of each outer ring instead
; of the last pixels of the inner ring for the boundary. That way, the inside of each 
; outer ring and the outside of each inner ring are ignored.
anglemax_arr = [360., 40., 35.] ; max. parallactic angle between ref. and science frame

corr_thresh_arr = [0.997, 0.9991]

obj='Alcor'
do_create_basis=1; need a new basis for each frame.

output_path = '/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/Alcor/macbook_' +$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)

; start to loop through parameters, where saving a new basis and binned cubes for
; for each new sz bin and/or filter size

foreach corr_thresh, corr_thresh_arr do begin

alcor_pipeline, pre_inj=1, neg_inj=0, coadd=coadd,$
								uncert=0, klip=0,use_gauss=0, fs=0, nod='dx_only',$
								cds='endcds', pre_clean=0, k_klip=k_klip, adi=0,$
								do_annmode=1, szbin=szbin, angsep=angsep,$
								bin_type=bin_type, filter=filter, n_ang=n_ang,$
								pa_off=pa_off,combine_type=combine_type,do_rdi=0,$
								annmode_inout=annmode_inout,do_dewarp=do_dewarp,$
								save_basis=save_basis, anglemax=anglemax,$
								rdi_extra=extra, extra=extra,$
								do_new_binning=do_new_binning,$
								do_create_basis=do_create_basis, lp_filter=lp_filter,$
								corr_thresh=corr_thresh, dxklip=0
								
ct = corr_thresh
; Alcor cubes
refcube1 = readfits(output_path+'/processed_right/dith1/'+obj+string(ct)+'_cube_skysub_cen_clean.fits')
refcube2 = readfits(output_path+'/processed_right/dith2/'+obj+string(ct)+'_cube_skysub_cen_clean.fits')
just_alcor_refcube = REFORM([REFORM(refcube1, (size(refcube1))[-1]), REFORM(refcube2, (size(refcube2))[-1])],$
	380,380,(size(refcube1))[3] + (size(refcube2))[3]);

; restore the Alcor cube angles (both nodding positions)
restore,filename=output_path + '/processed_right/dith1/' + obj + string(ct) +  '_parang_clean.sav'
just_alcor_ref_angles = angles

restore,filename=output_path + '/processed_right/dith2/' + obj + string(ct) +  '_parang_clean.sav'
just_alcor_ref_angles = [just_alcor_ref_angles, angles]

do_new_binning = 1
foreach szbin, bin_arr do begin

	do_new_binning = 1
	
	foreach bin_type, bin_type_arr do begin
		do_new_binning=1
		
		foreach combine_type, combine_type_arr do begin
			foreach n_ang, n_ang_arr do begin
				foreach filter, filter_arr do begin
				
					do_new_binning=1; already computed.
					
					foreach angsep, angsep_arr do begin
					
						foreach anglemax, anglemax_arr do begin
						
						foreach k_klip, kklip_arr do begin
						
							print, newline, 'running Alcor pipeline in loop for:',newline
							print, 'szbin = ', string(szbin)
							print, 'n_ang = ', string(n_ang)
							print, 'angsep = ', string(angsep)
							print, 'filter= ', string(filter)
							print, 'bin_type: ', bin_type
							print, 'combine_type:', combine_type
							print, 'k_klip = ', string(k_klip), newline
							print, 'do_new_binning: ', do_new_binning
							print, 'do_create_basis: ', do_create_basis, newline
							
							alcor_pipeline, pre_inj=0, neg_inj=0, coadd=coadd,$
								uncert=0, klip=0,use_gauss=0, fs=0, nod='dx_only',$
								cds='endcds', pre_clean=0, k_klip=k_klip, adi=0,$
								do_annmode=1, szbin=szbin, angsep=angsep,$
								bin_type=bin_type, filter=filter, n_ang=n_ang,$
								pa_off=pa_off,combine_type=combine_type,do_rdi=1,$
								annmode_inout=annmode_inout,do_dewarp=do_dewarp,$
								save_basis=save_basis, anglemax=anglemax,$
								rdi_extra=extra, extra=extra,$
								do_new_binning=do_new_binning,$
								do_create_basis=do_create_basis, lp_filter=lp_filter,$
								corr_thresh=corr_thresh, dxklip=0, refcube=just_alcor_refcube,$
								ref_ang=just_alcor_ref_angles
									
							do_new_binning=0; already made the basis for this loop (only filter, bin, will change it)
								
						endforeach; k_klip loop
						endforeach; anglemax loop
					endforeach; angsep loop
				endforeach; filter loop
			endforeach; nangarr
		endforeach; combine_type
	endforeach; bin_type loop
endforeach;szbin loop
endforeach; corr_thresh loop

end