pro kklip_loop_rdi5, kklip_arr=kklip_arr, angsep_arr=angsep_arr, bin_arr=bin_arr
;alcor_pipeline, pre_inj=1,pre_clean=1,coadd=25,klip=1,adi=0,fs=0,nod='dx_only',cds='endcds',k_klip=5,do_annmode=1,szbin=4,angsep=1,bin_type='median',filter=25.,n_ang=1,pa_off=0,combine_type='median'

;+
; NAME:
;       KKLIP_LOOP_RDI

; PURPOSE:
;       run alcor_pipeline.pro (or another reduction pipeline) in a loop with rdi=1
;		in order to test different RDI setups (e.g., with different k_klip or 
;		other kwarg values)

; CALLING SEQUENCE:
;       kklip_loop_rdi
;
; MODIFICATION HISTORY:
;       WRITTEN, March 2025 Gabriel Weible
;
; NOTE: if the reference cube used will be identical for all target frames
; (e.g., a true reference with no restriction on parallactic angle rotation
; like with cADI reference frames included), then immediately after the first target
; frame induces the creation of a basis (saved to the desktop by default) then the
; loop can be stopped and create_basis can be set to zero. Note that the basis
; even for an identical reference cube will also depend on bin and filter, so
; create_basis may be modified accordingly between changes in these values.
; Really, I should update rdiklip2 such that there is an option to create the basis
; for the first target frame, and *then* no longer create a new basis. This would
; be more dynamic for this sort of loop while saving on computation time.
;-

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)

new_binning_i=1

fwhm=13.0773; M-band
;lambda_over_D = 11.07 ; px
coadd=20
pa_off=0.
;save_basis=1; start by saving (only necessary as filter, bin change)
save_basis=0; don't save the basis for "RDI" that's really klip...
do_dewarp=1; don't dewarp if testing RDI

; this should roughly go fastest to slowest, except for the
; loop over kklip (which I think is important)
; I probably should loop over angsep deeper into the nested for loops
; so that this relatively important parameter may be for quickly varied, as well

kklip_arr=[6, 7, 8, 9, 10, 11, 12, 13, 14];[4,5,6,7];[4,5,6,7];new array from initial results;[8,9]; partially complete;[2, 3, 4, 5, 6, 7, 8, 9]
angsep_arr=[0.250, 0.255, 0.260, 0.265, 0.270, 0.275, 0.280, 0.285, 0.290, 0.295, 0.300, 0.305];[0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25];[0.1,0.2,0.3,0.4,0.5,0.8,1.6];[0.05, 0.06, 0.07, 0.7, 0.8, 0.9];[0.25/fwhm, 0.5/fwhm, 0.75/fwhm, 1.0/fwhm, 2.0/fwhm, 3.0/fwhm, 3.25/fwhm, 3.5/fwhm, 3.75/fwhm, 4.0/fwhm, 5.0/fwhm]

bin_type_arr = ['median']
combine_type_arr = ['nwadi']
filter_arr = [25.];[30., 40., 50., 60., 70., 80., 90.]
n_ang_arr = [2]
bin_arr = [5];I should probably try lower coadds as well!
annmode_inout=[8,34]; nwadi needs 1 extra px on either side?
anglemax_arr = [32., 34., 36., 38., 40., 42., 44., 46.] ; max. parallactic angle between ref. and science frame

corr_thresh = 0.92
ct = corr_thresh
obj='Alcor'
create_basis=1; need a new basis for each frame.

xp = 179.0-1.; back to no pa_off;189.163; approx. x-pixel position of Alcor B (w/pa_off=22.)
yp = 209.4-1.; back to no pa_off;198.093 ; approx. y-pixel position of Alcor B (w/pa_off=22.)

output_path = '/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/Alcor/macbook_' +$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)

save, filename=output_path+'/kklip_loop_arrays_rdi5.sav',kklip_arr, angsep_arr, filter_arr,$
	bin_type_arr, combine_type_arr, n_ang_arr, bin_arr

; Alcor cubes
refcube1 = readfits(output_path+'/processed_right/dith1/Alcor     0.920000_cube_skysub_cen_clean.fits')
refcube2 = readfits(output_path+'/processed_right/dith2/Alcor     0.920000_cube_skysub_cen_clean.fits')
just_alcor_refcube = REFORM([REFORM(refcube1, (size(refcube1))[-1]), REFORM(refcube2, (size(refcube2))[-1])],$
	380,380,(size(refcube1))[3] + (size(refcube2))[3]);

; restore the Alcor cube angles (both nodding positions)
restore,filename=output_path + '/processed_right/dith1/' + obj + string(ct) +  '_parang_clean.sav'
just_alcor_ref_angles = angles

restore,filename=output_path + '/processed_right/dith2/' + obj + string(ct) +  '_parang_clean.sav'
just_alcor_ref_angles = [just_alcor_ref_angles, angles]

; start to loop through parameters, where saving a new basis and binned cubes for
; for each new sz bin and/or filter size
foreach szbin, bin_arr do begin

	new_binning_i=1; already computed.
	;save_basis=1
	
	foreach bin_type, bin_type_arr do begin
		foreach combine_type, combine_type_arr do begin
			foreach n_ang, n_ang_arr do begin
				foreach filter, filter_arr do begin
				
					new_binning_i=1; already computed.
					;save_basis=1
					
					foreach angsep, angsep_arr do begin
					
						; I need these here for 'RDI' using science frames where
						; angsep actually needs to change.
						;save_basis=1
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
							
							alcor_pipeline, pre_inj=0, neg_inj=0, coadd=coadd, uncert=0, klip=0,$
							use_gauss=0, fs=0, nod='dx_only', cds='endcds', pre_clean=0,$
							k_klip=k_klip, adi=0, do_annmode=1, szbin=szbin, angsep=angsep,$
							bin_type=bin_type, filter=filter, n_ang=n_ang, pa_off=pa_off,$
							combine_type=combine_type,do_rdi=1,refcube=just_alcor_refcube,$
							annmode_inout=annmode_inout,create_basis=create_basis,$
							ref_ang=just_alcor_ref_angles,new_binning=new_binning_i,$
							do_dewarp=do_dewarp, save_basis=save_basis, anglemax=anglemax,$
							rdi_extra='5', extra='rdi5', mask_pt=[xp, yp]
								
							new_binning_i=0; already made the basis for this loop (only filter, bin, will change it)
							;save_basis=0
								
						endforeach; k_klip loop
						endforeach; anglemax loop
					endforeach; angsep loop
				endforeach; filter loop
			endforeach; nangarr
		endforeach; combine_type
	endforeach; bin_type loop
endforeach;szbin loop
end