pro snr_map, coadd=coadd, kklip_arr=kklip_arr, angsep_arr=angsep_arr,$
  bin_arr=bin_arr, extra=extra, filename=filename, save_map=save_map,$
  type=type, csv=csv, search=search
;+
; NAME:
;       KKLIP_LOOP_SNR

; PURPOSE:
;       Calculate SNR for reductions performed over various parameters (e.g., K_KLIP)

;       WRITTEN, 2025 Gabriel Weible
;       OPTIMIZED, 2025
;       PERFORMANCE OPTIMIZED, 2025
;-

compile_opt IDL2
t_start = systime(/seconds) ; For performance timing
newline = string(10B)

if not keyword_set(csv) then csv=0
if not keyword_set(save_map) then save_map=0

; USER INPUT
;-----------------------------------------------------------
fwhm = 13.382; M-band LMIRCam FWHM with dgvAPP360 (approx.)
boxehs = fwhm/2. ; source exclusion half-size
;nantr = 1000 ; stdev limit for masking remaining bad pixels
aper_rad = boxehs ; aperture radius for SNR calculation

; minimum and maximum radii for which to create the SNR map
; 7 seems like the innermost meaningful radius, and 38 is about the max
; for a total radius of 45 pixels with fwhm/2 left on the outside for the
; apertures.
r_min = 0 ; px
r_max = 170 ; px
dr = 0.25 ; px separation between noise/mean rings

;; Location of Alcor B to mask
;; (vip 91x91 sizing)
xp_mask = 166
yp_mask = 191

; how many azimuthal angles to calculate noise and signal at for each radius?
n_azimuth = 48

; Image size, HWHM Gaussian PSF for low-pass filtering
xs = 351 & ys = 351

; low- or high-pass filter
lp_width = 0
hp_width = 0

; number of sub-dr dithers in either dimension for each aperture
n_dith_1d = 9
;-----------------------------------------------------------

; Define center points and prepare arrays
xhs = xs/2.-0.5 & yhs = ys/2.-0.5 ; half-sizes

;; Pre-compute radius and angle of source to be masked (Alcor B)
rho_mask = sqrt((xhs-xp_mask)^2. + (yhs-yp_mask)^2.)
theta_mask = acos((xp_mask-xhs) / rho_mask)

output_path = '/Users/gweible/Library/CloudStorage/OneDrive - University of Arizona/research/Alcor/macbook_' +$
  strcompress(coadd, /r)
if keyword_set(extra) then output_path += '_'+strcompress(extra, /r)

if (type eq 'vip_klip') or (type eq 'vip_annklip') then output_path = '/Users/gweible/OneDrive - University of Arizona/research/Alcor/VIP'; search wherever the VIP files are.
if (type eq 'desktop') then output_path = '/Users/gweible/Desktop'

; File pattern matching - fix for the array error
if search eq 1 then begin

	if type eq 'klip' then begin
  		patterns = ['Alcor_k_klip_[0-9][0-9].0_comb_type_*_combined.fits', $ ; 2-digit
             'Alcor_k_klip_[0-9].00_comb_type_*_combined.fits', $      ; 1-digit
             'Alcor_k_klip_[0-9][0-9][0-9]._comb_type_*_combined.fits'] ; 3-digit
             
	endif
	;Alcor_keepnumber_700_filt_0.00000_neg_inj_0_uncert_0_szbin_1_type_mean_comb_type_median_normal_0_peak_thresh_0.960000_stddev_thresh_0.900000_right_adi.fits
	if type eq 'adi' then begin
		;patterns = ;[;'Alcor_ct_0.[0-9][0-9][0-9][0-9][0-9][0-9]_filt_0.00000_*_normal_0_right_adi.fits',$
		patterns =	['Alcor_keepnumber_[0-9][0-9][0-9]_filt_0.00000_*_normal_0_*_right_adi.fits',$
						'Alcor_keepnumber_[0-9][0-9][0-9][0-9]_filt_0.00000_*_normal_0_*_right_adi.fits']
	
		;patterns = ['Alcor_ct_0.[0-9][0-9][0-9][0-9][0-9][0-9]_filt_*_normal_0_right_adi.fits',$
		;			'Alcor_ct_0.[0-9][0-9][0-9][0-9][0-9][0-9]_filt_*_normal_1_right_adi.fits']
	endif
	
	if type eq 'vip_klip' then begin
		
		patterns =	['alcor_vip_klip_pcat_*_median_comb.fits']
		
	endif
	
	if type eq 'vip_annklip' then begin
	
		patterns = ['alcor_vip_klip_pcat_drot_*_annular_med_to_mean.fits',$
					'alcor_vip_klip_pcat_drot_*_annular_med_to_mean_2seg.fits']
		
			;patterns =	['ncomp_[0-9]_drot_0.[0-9]_mean_annular_frame.fits',$
			;				'ncomp_[0-9]_drot_0.[0-9][0-9]_mean_annular_frame.fits',$
			;				'ncomp_[0-9]_drot_0.[0-9][0-9][0-9]_mean_annular_frame.fits']
		
	endif
	
	if type eq 'desktop' then begin
		patterns = ['ncomp_*.fits']
	endif
  
  all_files = []
  for i=0, n_elements(patterns)-1 do begin
    files = FILE_SEARCH(output_path, '/'+patterns[i], COUNT=filecount)
    if filecount gt 0 then begin
      print, newline, 'Found ', filecount, ' files matching pattern ', i+1
      all_files = [all_files, files]
    endif
  endfor
endif else begin
  all_files = [filename]
endelse

; Pre-compute array of radii 
r_arr = findgen(fix(1./dr)*(r_max - r_min))*dr + r_min

; Pre-calculate dithering offsets (using mesh grid for efficiency)
xp_dith_offsets = (findgen(n_dith_1d) - (n_dith_1d-1)/2.) * (dr/2.)
yp_dith_offsets = xp_dith_offsets
n_dith = n_dith_1d * n_dith_1d

; Pre-compute ntheta values for each radius
ntheta_arr = fix(2.*!PI*r_arr/fwhm) > 3 ; Ensure at least 3

; Pre-compute stitheta and stendtheta for each radius (vectorized)
stitheta_arr = intarr(n_elements(r_arr))
for i = 0, n_elements(r_arr)-1 do begin
  stitheta_arr[i] = 1;ntheta_arr[i] lt 6 ? 1 : 2
endfor
stendtheta_arr = intarr(n_elements(r_arr))
n_theta_test_arr = intarr(n_elements(r_arr))

for i = 0, n_elements(r_arr)-1 do begin
  stendtheta_arr[i] = ntheta_arr[i] lt 6 ? ntheta_arr[i]-1 : ntheta_arr[i]-2
  n_theta_test_arr[i] = (stendtheta_arr[i] - stitheta_arr[i] + 1) > 1
endfor

; Pre-compute all azimuthal angles
az_angles = findgen(n_azimuth) * 2*!PI / n_azimuth

; Pre-compute coordinate grid once
xx_grid = rebin(indgen(fix(xs)), fix(xs), fix(ys))
yy_grid = rebin(reform(indgen(fix(ys)), 1, fix(ys)), fix(xs), fix(ys))
rho_grid = sqrt((xx_grid - xhs)^2 + (yy_grid - yhs)^2)

; Prepare array for peak results
peak_results = fltarr(n_elements(all_files), /nozero)

for filei = 0, n_elements(all_files)-1 do begin
    file = all_files[filei]
    print, 'Processing file:', newline, file
    
    ; Read and process image (just once)
    zimage = readfits_fast(file)
    
    ; High-pass filter
    if hp_width gt 0 then zimage -= filter_image(zimage, fwhm=hp_width,$
    	psf=psf_hp, /all_pixels)
    ; Low-pass filter
    if lp_width gt 0 then zimage = filter_image(zimage, fwhm=lp_width,$
    	psf=psf_lp, /all_pixels)
    
    ; Create SNR image with same dimensions
    SNR_img = fltarr(fix(xs), fix(ys), /nozero)

	; NOISE IMAGE SETUP
    nimage = zimage
    
  ;  ; Mask Alcor B in noise image - vectorized
  ;  mask_indices = where(sqrt((xx_grid - xp_mask)^2 + (yy_grid - yp_mask)^2) le boxehs, mask_count)
  ;  nimage[mask_indices] = !values.f_nan
    
     ; Convert NaNs to zeros for aperture photometry
    nan_indices = where(finite(zimage) ne 1, nan_count)
    if nan_count gt 0 then zimage[nan_indices] = median(zimage, /double, /even)
    
    nan_indices = where(finite(nimage) ne 1, nan_count)
    if nan_count gt 0 then nimage[nan_indices] = median(nimage, /double, /even)
    
    if keyword_set(save_map) then begin
    	print, 'saving filtered image...'
    	print, 'output_path: ', output_path
    	print, 'file: ', file
    	writefits, output_path + '/SNR_maps/' + file_basename(file, '.fits') + '_filtered.fits', zimage
    
	endif
	
    ; Pre-allocate arrays for noise and mean
    noise_rhos = fltarr(n_elements(r_arr), /nozero)
    mean_rhos = fltarr(n_elements(r_arr), /nozero)
    
    ; For better performance, pre-calculate all test positions
    all_test_x = fltarr(n_elements(r_arr), n_azimuth, n_dith, max(n_theta_test_arr),$
    	/nozero)
    all_test_y = fltarr(n_elements(r_arr), n_azimuth, n_dith, max(n_theta_test_arr),$
    	/nozero)
    all_valid = bytarr(n_elements(r_arr), n_azimuth, n_dith, max(n_theta_test_arr))
    
    ; Calculate all test positions in advance
    for r_idx = 0, n_elements(r_arr)-1 do begin
        rho_px = r_arr[r_idx]
        ;print, 'rho_px:', rho_px
        ntheta = ntheta_arr[r_idx]
        stitheta = stitheta_arr[r_idx]
        n_theta_test = n_theta_test_arr[r_idx]
        
        for az_idx = 0, n_azimuth-1 do begin
            dtheta_az = az_angles[az_idx]
            
            ; Base position
            xp_base = rho_px * cos(theta_mask + dtheta_az) + xhs
            yp_base = rho_px * sin(theta_mask + dtheta_az) + yhs
            
            dith_idx = 0
            for x_off_idx = 0, n_dith_1d-1 do begin
                for y_off_idx = 0, n_dith_1d-1 do begin
                    xp_dith = xp_base + xp_dith_offsets[x_off_idx]
                    yp_dith = yp_base + yp_dith_offsets[y_off_idx]
                    
                    rho_dith = sqrt((xhs-xp_dith)^2 + (yhs-yp_dith)^2)
                    theta_dith = acos((xp_dith-xhs) / rho_dith)
                    
                    for itheta_idx = 0, n_theta_test-1 do begin
                        itheta_dith = stitheta + itheta_idx
                        ttheta_dith = theta_dith + (float(itheta_dith)/float(ntheta))*2.*!PI
                        
                        xtest = xhs + (float(rho_dith)*cos(ttheta_dith))
                        ytest = yhs + (float(rho_dith)*sin(ttheta_dith))
                        
                        ; Check bounds
                        is_valid = (xtest ge 0) and (xtest lt xs) and (ytest ge 0) and (ytest lt ys)
                        
                        all_test_x[r_idx, az_idx, dith_idx, itheta_idx] = xtest
                        all_test_y[r_idx, az_idx, dith_idx, itheta_idx] = ytest
                        all_valid[r_idx, az_idx, dith_idx, itheta_idx] = is_valid
                    endfor
                    
                    dith_idx++
                endfor
            endfor
        endfor
    endfor
    
    ; Process each radius more efficiently by grouping aperture photometry
    for r_idx = 0, n_elements(r_arr)-1 do begin
        n_theta_test = n_theta_test_arr[r_idx]
        flux_measurements = fltarr(n_azimuth, n_dith, n_theta_test, /nozero)
        
        ; Group all aperture measurements for this radius
        for az_idx = 0, n_azimuth-1 do begin
            for dith_idx = 0, n_dith-1 do begin
                for itheta_idx = 0, n_theta_test-1 do begin
                    if all_valid[r_idx, az_idx, dith_idx, itheta_idx] then begin
                        xtest = all_test_x[r_idx, az_idx, dith_idx, itheta_idx]
                        ytest = all_test_y[r_idx, az_idx, dith_idx, itheta_idx]
                        
                        ; Use existing aperture function, but could be optimized further
                        aper, nimage, xtest, ytest, flux, fluxerr, sky, $
                            skyerr, 1, aper_rad, [0,0], [-99E99,99E99], /flux, /silent, SETSKYVAL=0
                        
                        flux_measurements[az_idx, dith_idx, itheta_idx] = flux
                    endif else begin
                        flux_measurements[az_idx, dith_idx, itheta_idx] = !values.f_nan
                    endelse
                endfor
            endfor
        endfor
        
        ; Zero fluxes to NaN
        zero_indices = where(flux_measurements eq 0., zero_count)
        if zero_count gt 0 then flux_measurements[zero_indices] = !values.f_nan
        
        ; Process measurements for this radius
        az_noise = fltarr(n_azimuth, /nozero)
        az_mean = fltarr(n_azimuth, /nozero)
        
        for az_idx = 0, n_azimuth-1 do begin
            dith_noise = fltarr(n_dith, /nozero)
            dith_mean = fltarr(n_dith, /nozero)
            
            for dith_idx = 0, n_dith-1 do begin
                fluxes = flux_measurements[az_idx, dith_idx, *]
                valid_flux = where(finite(fluxes), n_valid_flux)
                
                if n_valid_flux gt 0 then begin
                    dith_noise[dith_idx] = stddev(fluxes[valid_flux], /double) * sqrt(1.+(1./n_valid_flux))
                    dith_mean[dith_idx] = mean(fluxes[valid_flux], /double)
                endif else begin
                    dith_noise[dith_idx] = !values.f_nan
                    dith_mean[dith_idx] = !values.f_nan
                endelse
            endfor
            
            valid_dith = where(finite(dith_noise), n_valid_dith)
            if n_valid_dith gt 0 then begin
                az_noise[az_idx] = mean(dith_noise[valid_dith], /double)
                az_mean[az_idx] = mean(dith_mean[valid_dith], /double)
            endif else begin
                az_noise[az_idx] = !values.f_nan
                az_mean[az_idx] = !values.f_nan
            endelse
        endfor
        
        valid_az = where(finite(az_noise), n_valid_az)
        if n_valid_az gt 0 then begin
            noise_rhos[r_idx] = mean(az_noise[valid_az], /double)
            mean_rhos[r_idx] = mean(az_mean[valid_az], /double)
        endif else begin
            noise_rhos[r_idx] = !values.f_nan
            mean_rhos[r_idx] = !values.f_nan
        endelse
    endfor
    
    ; Find valid pixels for SNR calculation - vectorized approach
    valid_mask = (rho_grid ge r_min) and (rho_grid le r_max) and finite(zimage)
    valid_indices = where(valid_mask, count_valid)
    
    if count_valid gt 0 then begin
        ; Extract coordinates
        valid_x = xx_grid[valid_indices]
        valid_y = yy_grid[valid_indices]
        valid_rho = rho_grid[valid_indices]
        
        ; Clean noise and mean arrays for interpolation
        valid_r_arr = where(finite(noise_rhos), n_valid_noise)
        if n_valid_noise gt 1 then begin
            r_arr_clean = r_arr[valid_r_arr]
            noise_rhos_clean = noise_rhos[valid_r_arr]
            mean_rhos_clean = mean_rhos[valid_r_arr]
            
            ; Process in larger chunks for better performance
            chunk_size = 50000L ; Larger chunks for fewer iterations
            for chunk_start = 0L, count_valid-1, chunk_size do begin
                chunk_end = (chunk_start + chunk_size - 1) < (count_valid - 1)
                chunk_indices = lindgen(chunk_end - chunk_start + 1) + chunk_start
                
                ; Get chunk coordinates
                chunk_x = valid_x[chunk_indices]
                chunk_y = valid_y[chunk_indices]
                chunk_rho = valid_rho[chunk_indices]
                
                ; Interpolate noise and mean for entire chunk
                chunk_noise = INTERPOL(noise_rhos_clean, r_arr_clean, chunk_rho, /spline)
                chunk_mean = INTERPOL(mean_rhos_clean, r_arr_clean, chunk_rho, /spline)
                
                ; Batch aperture photometry if possible
                chunk_signal = fltarr(n_elements(chunk_indices), /nozero)
                
                ; Unfortunately we still need to do aperture photometry one by one
                ; This could be further optimized with a custom aperture function
                for i = 0L, n_elements(chunk_indices)-1 do begin
                    aper, zimage, chunk_x[i], chunk_y[i], signal, err, sky, skyerr, 1, $
                        aper_rad, [0,0], [-99E99,99E99], /flux, /silent, SETSKYVAL=0
                    
                    chunk_signal[i] = signal
                endfor
                
                ; Calculate SNR for entire chunk
                chunk_snr = (chunk_signal - chunk_mean) / chunk_noise
                
                ; Store results in SNR image
                SNR_img[valid_indices[chunk_indices]] = chunk_snr
            endfor
        endif
    endif
    
    ; Fix NaNs and apply final smoothing
    nan_indices = where(finite(SNR_img) ne 1, nan_count)
    if nan_count gt 0 then SNR_img[nan_indices] = 0.
    
    for xx=0,xs-1 do begin
    for yy=0,ys-1 do begin
    	if sqrt( (( xx - xhs)^2.) + (( yy - yhs)^2.) ) le r_min then SNR_img[xx,yy,*] = 0.
		if sqrt( (( xx - xhs)^2.) + (( yy - yhs)^2.) ) ge r_max then SNR_img[xx,yy,*] = 0.
    endfor
    endfor
    
    ; Low-pass filter
    if lp_width gt 0 then SNR_img = filter_image(SNR_img, fwhm=lp_width,$
    	psf=psf_lp, /all_pixels)
    
    ; Find peak SNR
    peak_SNR = max(SNR_img, max_loc)
    print, newline, 'peak_SNR:', peak_SNR, newline
    peak_results[filei] = peak_SNR
    
    ; Optional: write individual SNR maps if needed
    if save_map eq 1 then begin
    	print, 'saving SNR map image...'
    	print, 'output_path: ', output_path
    	print, 'file: ', file
    	writefits, output_path + '/SNR_maps/' + file_basename(file, '.fits') + '_SNR_map.fits', SNR_img
    endif
endfor

; Write results to CSV
if csv eq 1 then WRITE_CSV, output_path + '/Alcor_contrast_output_big_' + type + '.csv', $
    peak_results, all_files

t_end = systime(/seconds)
print, 'Processing complete for all files.'
print, 'Total execution time: ', t_end - t_start, ' seconds'
end