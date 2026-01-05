pro tyc5709_klip_loop

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)

fwhm=10.6245 ; from VIP (L band)
coadd=1; # really, I'll use data pre-co-added by VIP.
pa_off=0.
;save_basis=1; start by saving (only necessary as filter, bin change)
save_basis=0; don't save the basis for "RDI" that's really klip...
do_dewarp=1

; this should roughly go fastest to slowest, except for the
; loop over kklip (which I think is important)
; I probably should loop over angsep deeper into the nested for loops
; so that this relatively important parameter may be for quickly varied, as well
kklip_arr = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
angsep_arr=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]

; only keep_number cross-correlation cleaning...
stddev_thresh=0
peak_thresh=0

bin_type_arr = ['res_mean']
combine_type_arr = ['median']
filter_arr = [0.0]
n_ang_arr = [1]
szbin = 1

; this should be 2x FWHM in width after everything is actually combined
annmode_inout=[0, 60];
annmode_in = annmode_inout[0]
annmode_out = annmode_inout[1]
r_px = ((annmode_out + annmode_in) / 2.)

; 1-px of overlap. Then, I just take 2nd-from-inside-most of each outer ring instead
; of the last pixels of the inner ring for the boundary. That way, the inside of each 
; outer ring and the outside of each inner ring are ignored.
anglemax_arr = [360.] ; max. parallactic angle between ref. and science frame

obj='tyc5709'
do_create_basis=1; need a new basis for each frame.

output_path = '/Volumes/T7/TYC5709/gabe_macbook' +$
	strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)

; start to loop through parameters, where saving a new basis and binned cubes for
; for each new sz bin and/or filter size
do_new_binning = 1
foreach n_ang, n_ang_arr do begin
    foreach filter, filter_arr do begin
    
        do_new_binning=1; already computed.
        
        foreach angsep, angsep_arr do begin
            do_new_binning=1
        
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
                    pre_clean=0, k_klip=k_klip, adi=0,$
                    do_annmode=1, szbin=szbin, angsep=angsep,$
                    bin_type=bin_type, filter=filter, n_ang=n_ang,$
                    pa_off=pa_off,combine_type=combine_type,do_rdi=1,$
                    annmode_inout=annmode_inout,do_dewarp=do_dewarp,$
                    save_basis=save_basis, anglemax=anglemax,$
                    rdi_extra=extra, extra=extra,$
                    do_new_binning=do_new_binning,$
                    do_create_basis=do_create_basis,$
                    dxklip=1,$
                    stddev_thresh=stddev_thresh, peak_thresh=peak_thresh
                    
                tyc5709_pipeline, pre_inj=0, neg_inj=0, coadd=coadd,$
                    use_gauss=0, klip=0, fs=0, nod=5, pre_clean=0,$
                    k_klip=k_klip, adi=0, do_annmode=1, szbin=szbin,$
                    angsep=angsep, bin_type=bin_type, filter=filter,$
                    pa_off=pa_off, n_ang=n_ang, combine_type=combine_type,$
                    do_rdi=1, annmode_inout=annmode_inout, do_new_binning=do_new_binning,$
                    do_create_basis=do_create_basis, do_dewarp=do_dewarp,$
                    save_basis=save_basis, anglemax=anglemax, rdi_extra=extra,$
                    dxklip=1, peak_thresh=peak_thresh, stddev_thresh=stddev_thresh
                        
                do_new_binning=0; already made the basis for this loop (only filter, bin, will change it)
                    
            endforeach; k_klip loop
            endforeach; anglemax loop
        endforeach; angsep loop
    endforeach; filter loop
endforeach; nangarr

end
