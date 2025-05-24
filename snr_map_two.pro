;+
; NAME:
;   snr_map_two
; PURPOSE:
;   Compute an SNR map from one or two images by pooling noise apertures and combining signals.
;+
pro snr_map_two, filename1, filename2, fwhm=fwhm, r_min=r_min, r_max=r_max, dr=dr, $
    n_az=n_az, boxehs=boxehs, lp_width=lp_width, hp_width=hp_width, save_map=save_map
    
    compile_opt IDL2

    ;-- set keyword defaults
    if not keyword_set(boxehs) then boxehs = 6.53865
    if not keyword_set(fwhm)   then fwhm   = 13.0773
    if not keyword_set(r_min)  then r_min  = 0.
    if not keyword_set(r_max)  then r_max  = 180.
    if not keyword_set(dr)     then dr     = 0.1
    if not keyword_set(n_az)   then n_az   = 36
    if not keyword_set(lp_width) then lp_width = 0.
    if not keyword_set(hp_width) then hp_width = 0.
    if not keyword_set(save_map) then save_map = 1

    ;-- detect if second image provided
    do_two = (strlen(filename2) gt 0)

    ;-- aperture radius
    aper_rad = boxehs

    ;-- read images
    img1 = readfits_fast(filename1)
    if do_two then img2 = readfits_fast(filename2)

    ;-- replace NaNs with median
    idx = where(finite(img1) eq 0, n)
    if n gt 0 then img1[idx] = median(img1, /even)
    if do_two then begin
        idx = where(finite(img2) eq 0, n)
        if n gt 0 then img2[idx] = median(img2, /even)
    endif

    ;-- high-pass filter
    if hp_width gt 0 then begin
        hp_npix = floor(7 * hp_width)
        if ~ODD(hp_npix) then hp_npix++
        psf_hp = psf_Gaussian(NPIX=hp_npix, FWHM=[hp_width, hp_width], /normalize)
        img1 -= convolve(img1, psf_hp)
        if do_two then img2 -= convolve(img2, psf_hp)
    endif

    ;-- low-pass filter
    if lp_width gt 0 then begin
        lp_npix = floor(7 * lp_width)
        if ~ODD(lp_npix) then lp_npix++
        psf_lp = psf_Gaussian(NPIX=lp_npix, FWHM=[lp_width, lp_width], /normalize)
        img1 = convolve(temporary(img1), psf_lp)
        if do_two then img2 = convolve(temporary(img2), psf_lp)
    endif

    ;-- coordinate grid
    dims = size(img1, /dimensions)
    nx = dims[0]
    ny = dims[1]
    xhs = (nx-1)/2.
    yhs = (ny-1)/2.
    xx = rebin(indgen(nx), nx, ny)
    yy = rebin(reform(indgen(ny),1,ny), nx, ny)
    rho = sqrt((xx-xhs)^2 + (yy-yhs)^2)

    ;-- radii and noise array
    n_r = floor((r_max - r_min)/dr) + 1
    r_arr = findgen(n_r)*dr + r_min
    noise_rhos = fltarr(n_r)
    angles = findgen(n_az)*2*!PI/n_az

    ;-- compute noise by pooling
    for i=0,n_r-1 do begin
        count = n_az * (do_two ? 2 : 1)
        fluxs = fltarr(count)
        for j=0,n_az-1 do begin
            th = angles[j]
            x0 = xhs + r_arr[i]*cos(th)
            y0 = yhs + r_arr[i]*sin(th)
            aper, img1, x0, y0, f1, err, sk, skerr, 1, aper_rad, [0,0],$
                 [-1e10,1e10], /flux, /silent, /exact, SETSKYVAL=0
            fluxs[j] = f1
            if do_two then begin
                aper, img2, x0, y0, f2, err, sk2, sk2err, 1, aper_rad, [0,0],$
                     [-1e10,1e10], /flux, /silent, /exact, SETSKYVAL=0
                fluxs[j + n_az] = f2
            endif
        endfor
        noise_rhos[i] = stddev(fluxs)
    endfor

    ;-- build SNR map
    SNR = fltarr(nx, ny)
    valid = where((rho ge r_min) and (rho le r_max), nvalid)
    for k=0,nvalid-1 do begin
        idx = valid[k]
        x0 = xx[idx]
        y0 = yy[idx]
        r0 = rho[idx]
        ib = long((r0 - r_min)/dr + 0.5)
        ib = max([0, min([ib, n_r-1])])
        sigma = noise_rhos[ib]

        aper, img1, x0, y0, s1, err, sk, skerr, 1, aper_rad, [0,0],$
             [-1e10,1e10], /flux, /silent, /exact, SETSKYVAL=0
        if do_two then begin
            aper, img2, x0, y0, s2, err, sk2, sk2err, 1, aper_rad, [0,0],$
                 [-1e10,1e10], /flux, /silent, /exact, SETSKYVAL=0
            SNR[idx] = (s1 + s2) / (sigma * sqrt(2.))
        endif else begin
            SNR[idx] = s1 / sigma
        endelse
    endfor

    ;-- write if requested
    if save_map then begin
        out = file_basename(filename1,'.fits')
        if do_two then out += '_' + file_basename(filename2,'.fits')
        out += '_SNR_map.fits'
        writefits, out, SNR
        print, 'Wrote SNR map to ', out
    endif

    ;-- report peak
    peak = max(SNR, loc)
    print, 'Peak SNR = ', peak
end
