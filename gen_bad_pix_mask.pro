pro gen_bad_pix_mask, dark, flat, mask_out=mask_out, $
                               hot_sigma=hot_sigma, flat_thresh=flat_thresh,$
                               debug=debug;,$
                              ; hp_flat_thresh=hp_flat_thresh, hp_dark_thresh=hp_dark_thresh

    compile_opt idl2

    ; --- User-defined thresholds (optional keyword overrides) ---
    if n_elements(hot_sigma) eq 0 then hot_sigma = 3.0
    if n_elements(flat_thresh) eq 0 then flat_thresh = [0.48, 1.15]  ; [low, high]
    ; % deviation from local median behavior in re-normalized, high-pass filtered flat.
   ; if n_elements(hp_flat_thresh) eq 0 then hp_flat_thresh = [-0.1, 0.1]  ; [low, high]
    ; % deviation from local median behavior in re-normalized, high-pass filtered dark.
    ;if n_elements(hp_dark_thresh) eq 0 then hp_dark_thresh = [-500, 500]  ; [low, high]

	; RAW DARK:
    ; --- Statistics on dark frame ---
    ;dark_med = median(dark, /even, /double)
    ;dark_std = stddev(dark, /nan, /double)
    ; --- Identify bad pixels in dark frame ---
    ;hot_mask  = (dark GT (dark_med + hot_sigma * dark_std))
    ;dead_mask = (dark LT (dark_med - hot_sigma * dark_std))
    
    ; (PRE-NORMALIZED) FLAT:
     ; --- Identify bad pixels in flat field ---
     flat = flat / median(flat, /double, /even)
    low_resp_mask  = (flat LT flat_thresh[0])
    high_resp_mask = (flat GT flat_thresh[1])
    
    ; FILTERED FLAT (COMES NORMALIZED):
    ; high-pass filter flat field to identify local outliers
    flat_hp = flat - filter_image(flat, median=32, /all_pixels)
    ; re-normalize after HP filter
    norm_flat_hp = flat_hp / median(flat_hp, /even, /double)
    ; --- Identify bad pixels in the filtered + renormalized flat field ---
    ;low_resp_mask_hp  = (norm_flat_hp LT hp_flat_thresh[0])
    ;high_resp_mask_hp = (norm_flat_hp GT hp_flat_thresh[1])
    
    norm_flat_hp_med = median(norm_flat_hp, /even, /double)
    norm_flat_hp_std = stddev(norm_flat_hp, /nan, /double)
    ; --- Identify bad pixels in norm_flat_hp frame ---
    hot_norm_flat_hp_mask  = (norm_flat_hp GT (norm_flat_hp_med + hot_sigma * norm_flat_hp_std))
    dead_norm_flat_hp_mask = (norm_flat_hp LT (norm_flat_hp_med - hot_sigma * norm_flat_hp_std))
    
    ; FILTERED DARK
    ; normalize and high-pass filter dark field to identify local outliers
    dark_norm = dark / median(dark, /even, /double)
    dark_norm_hp = dark_norm - filter_image(dark_norm, median=32, /all_pixels)
    ; re-normalize after HP filter
    norm_dark_norm_hp = dark_norm_hp / median(dark_norm_hp, /even, /double)
    ; --- Identify bad pixels in filtered dark field ---
    ;hp_dead_mask  = (norm_dark_norm_hp LT hp_dark_thresh[0])
    ;hp_hot_mask  = (norm_dark_norm_hp GT hp_dark_thresh[1])
    
    norm_dark_norm_hp_med = median(norm_dark_norm_hp, /even, /double)
    norm_dark_norm_hp_std = stddev(norm_dark_norm_hp, /nan, /double)
    ; --- Identify bad pixels in norm_dark_norm_hp frame ---
    hot_norm_dark_norm_hp_mask  = (norm_dark_norm_hp GT (norm_dark_norm_hp_med + hot_sigma * norm_dark_norm_hp_std))
    dead_norm_dark_norm_hp_mask = (norm_dark_norm_hp LT (norm_dark_norm_hp_med - hot_sigma * norm_dark_norm_hp_std))

    ; --- Combine all masks ---
    ;bad_pixel_mask = hot_mask OR dead_mask OR low_resp_mask OR high_resp_mask OR $
    ;	low_resp_mask_hp OR high_resp_mask_hp OR hp_dead_mask OR hp_hot_mask
    
    bad_pixel_mask = low_resp_mask OR high_resp_mask OR hot_norm_flat_hp_mask OR $
    	dead_norm_flat_hp_mask OR hot_norm_dark_norm_hp_mask OR dead_norm_dark_norm_hp_mask

    ; --- Final good-pixel mask: 1 = good, 0 = bad ---
    good_pixel_mask = 1b - byte(bad_pixel_mask)
    
    if debug eq 1 then begin
		; test which of these is the problem...
		hot_norm_flat_hp_mask  = 1b - byte(hot_norm_flat_hp_mask)
		writefits, '~/Desktop/hot_norm_flat_hp_mask_test.fits', hot_norm_flat_hp_mask
		dead_norm_flat_hp_mask  = 1b - byte(dead_norm_flat_hp_mask)
		writefits, '~/Desktop/dead_norm_flat_hp_mask_test.fits', dead_norm_flat_hp_mask
		
		low_resp_mask  = 1b - byte(low_resp_mask)
		writefits, '~/Desktop/low_resp_mask_test.fits', low_resp_mask
		high_resp_mask  = 1b - byte(high_resp_mask)
		writefits, '~/Desktop/high_resp_mask_test.fits', high_resp_mask
		
		hot_norm_dark_norm_hp_mask  = 1b - byte(hot_norm_dark_norm_hp_mask)
		 writefits, '~/Desktop/hot_norm_dark_norm_hp_mask_test.fits', hot_norm_dark_norm_hp_mask
		dead_norm_dark_norm_hp_mask  = 1b - byte(dead_norm_dark_norm_hp_mask)
		writefits, '~/Desktop/dead_norm_dark_norm_hp_mask_test.fits', dead_norm_dark_norm_hp_mask
    
    endif

    ; Return output if requested
    if arg_present(mask_out) then mask_out = good_pixel_mask

end
