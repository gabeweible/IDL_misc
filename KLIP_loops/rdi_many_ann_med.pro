pro rdi_many_ann_med, kklip_arr=kklip_arr, angsep_arr=angsep_arr, bin_arr=bin_arr

compile_opt IDL2
newline = string(10B)

; ONLY NEEDS TO RUN AS = 1 ONCE IF FILTER, BINNING, DESTRIPING REMAIN THE SAME!
;-------------------------------------------------
do_new_binning = 1
;-------------------------------------------------

extra = 'rdi_many_ann_med'
fwhm=13.0773; M-band

image_sz_x = 380.
image_sz_y = 380.

strip_sz = 0

; star center
x_center = image_sz_x / 2.
y_center = image_sz_y / 2.

if x_center eq y_center then half_sz = x_center

; OPTIONS FOR OVERLAPPING ANNULI TO MEDIAN-COMBINE
;-------------------------------------------------
wr = 187;fix(7*fwhm);29;fix(2*fwhm); px annulus width. I'm trying 1x FWHM to speed things up.
inner_ann_start = 3;7; innermost radius considered (even numbers work better?)
ann_overlap = 1;6; just overlap 2 pixels?

; remaining radius to cover divided by shift of each annulus (wr - ann_overlap)
;num_annuli = fix((half_sz - inner_ann_start) / (wr - ann_overlap))
num_annuli = 1;4

; with the outermost bad pixel radius of the inner ring (this should realistically
; be the minimum overlap acceptable)
;-------------------------------------------------

; COMPUTE INNER AND OUTER ANNULI RADII
;-------------------------------------------------
annmode_ins = fltarr(num_annuli) + inner_ann_start ; start 8 annuli all off at 1st-annulus inner radius
annmode_ins_adds = findgen(num_annuli) * (wr - ann_overlap); shift each
; sequential annulus by the width of the rings, minus the desired overlap
annmode_ins += annmode_ins_adds; add the shifts to get all the sequential
; inner radii

; outer radii are just the inner radii + the ring width
annmode_outs = annmode_ins + wr

; loop through n_annuli to stitch together the inner and outer radii
annmode_inout_arr = [[]]
for j=0,num_annuli-1 do begin
	; append each [inner, outer] into a 2-D arr
	annmode_inout_arr = [[annmode_inout_arr], [annmode_ins[j], annmode_outs[j]] ]
endfor

;-------------------------------------------------

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)

;lambda_over_D = 11.07 ; px
coadd=20
pa_off=0.
;save_basis=1; start by saving (only necessary as filter, bin change)
save_basis=0; don't save the basis for "RDI" that's really klip...
do_dewarp=1; don't dewarp if testing RDI
neg_inj=0

pxscale_dx = 0.010610

; this should roughly go fastest to slowest, except for the
; loop over kklip (which I think is important)
; I probably should loop over angsep deeper into the nested for loops
; so that this relatively important parameter may be for quickly varied, as well

bin_type = 'median'
combine_type = 'median'
filter = 29.
n_ang = 2
szbin = 10
anglemax = 34.
nrings = 1

corr_thresh = 0.92
ct = corr_thresh
obj='Alcor'
do_create_basis=1; need a new basis for each frame.

output_path = '/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/Alcor/macbook_' +$
	strcompress(coadd, /r)
output_path += '_'+strcompress(extra, /r)
output_path += '/'

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

; loop through all the annuli desired to median-combine
big_3d_cube = list()
for j=0,num_annuli-1 do begin
							
	annmode_inout = annmode_inout_arr[*,j]; grab the desired annulus
	; adjust angsep with radius: arclength goes like radius, so if we want the actual
	; same references, then we should enforce a stricter angsep proportional to r
	
	; 21. is the center radius of our target annulus for Alcor B (annmode_inout_B = [8, 34])
	annmode_in = annmode_inout[0]
	annmode_out = annmode_inout[1]
	
	r_px = ((annmode_out + annmode_in) / 2.)
	r_arcsec = r_px * pxscale_dx
	angsep = 0.275 * r_px / 21.
	
	; trying variable k_klip
	k_klip=11
	;k_klip = fix( 7. + 0.1*EXP(4.-r_arcsec) + (0.075/r_arcsec^3.) - (0.31/r_arcsec^2.))
	
	; need for output filename to read it in for safekeeping before it's over-written in the loop
	super_suffix = strcompress(output_path + 'combined/' + obj +'_k_klip_'+$
		string(sigfig(k_klip,3)) +'_comb_type_'+combine_type+'_angsep_'+$
		string(sigfig(angsep,4))+'_angmax_'+string(sigfig(anglemax,3))+'_nrings_'+$
		string(sigfig(nrings, 2))+'_nang_'+string(sigfig(n_ang,2)) + '_filter_' +$
		string(filter) + '_bin_' + string(sigfig(szbin,2)) + '_type_' + bin_type +$
		'_neg_inj_' + string(neg_inj), /r)
	if keyword_set(trial) then super_suffix += '_trial_' + string(trial, format='(F0.4)')

	
	alcor_pipeline, pre_inj=0, neg_inj=neg_inj, coadd=coadd,$
		uncert=0, klip=0,use_gauss=0, fs=0, nod='dx_only',$
		cds='endcds', pre_clean=0, k_klip=k_klip, adi=0,$
		do_annmode=1, szbin=szbin, angsep=angsep,$
		bin_type=bin_type, filter=filter, n_ang=n_ang,$
		pa_off=pa_off,combine_type=combine_type,do_rdi=1,$
		refcube=just_alcor_refcube,$
		annmode_inout=annmode_inout,$
		ref_ang=just_alcor_ref_angles,do_dewarp=do_dewarp,$
		save_basis=save_basis, anglemax=anglemax,$
		rdi_extra='many_ann_med', extra=extra,$
		do_new_binning=do_new_binning, do_create_basis=do_create_basis
			
	do_new_binning=0; already made the basis for this loop (only filter, bin, will change it)
	
	; Write combined and individual nod results
	output_annulus = readfits(strcompress(super_suffix + '_combined.fits', /rem))
	
	; crop the first and last pixel radii which are noisy and weird.
	
	; Loop through all pixels in the image
	for x = 0, image_sz_x - 1 do begin
		for y = 0, image_sz_y - 1 do begin
		
			;Pythagorean theorem baby! (everything within the radius of pixels is turned into NaN)
			if sqrt( ((x-x_center)^2.) + ((y-y_center)^2.) ) le (annmode_in+1+strip_sz) then output_annulus[x,y] = !values.f_nan
			if sqrt( ((x-x_center)^2.) + ((y-y_center)^2.) ) ge (annmode_out-1-strip_sz) then output_annulus[x,y] = !values.f_nan
			
		endfor; y-pixel loop
	endfor; x-pixel loop
	
	big_3d_cube.Add, [[output_annulus]]; add 2-D annulus image to 3-D list
	
	; Writefits takes an array as input, so we'll need to convert our list over to an array
	print, 'Converting to array...'
	obj_cube = big_3d_cube.toArray(/TRANSPOSE)
	
	print, 'Writing cube FITS...'
	;Write the cube
	writefits, strcompress(super_suffix + '_3d_cube.fits',/rem), obj_cube
	print, 'Cube FITS written!'
	
	if j gt 0 then begin
		print, newline, 'Taking median...'
		running_medframe = median(obj_cube, dimension=3, /even, /double)
		running_mean = mean(obj_cube, dim=3, /double, /nan)
		
		print, 'Writing running medframe FITS...'
		;Write the cube
		writefits, strcompress(super_suffix + '_running_medframe.fits',/rem), running_medframe
		writefits, strcompress(super_suffix + '_running_mean.fits',/rem), running_mean
		print, 'Running medframe FITS written!', newline
	endif
	
endfor; loop over annuli indexed by j

end