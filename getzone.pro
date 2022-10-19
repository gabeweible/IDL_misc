function getzone, cube, annulus, angles

; AD
; Dec 2015 - Undisclosed Location
;
;
; The .indices returns an array of the (linear) indices of the pixels
; that are extracted in the zone. This can be used to copy these back,
; after processing, to the final image.
; =======================

 sc = size(cube)

 mask  = fltarr(sc[1],sc[2])

 dist_circle,  mask, [ sc[1],sc[2]], sc[1]/2.-1, sc[2]/2.-1

; xx = replicate(findgen(sc[1]),sc[2])-sc[1]/2.-1
; yy = transpose(replicate(findgen(sc[2]),sc[1]))-sc[2]/2.-1
 
 ang = mask
 ang = angmask(mask, sc[1]/2.-1, sc[2]/2.-1)

; good =mask
; good[*] = !values.f_nan
 good = where(ang ge angles[0]*!DTOR and ang le angles[1]*!DTOR and mask ge annulus[0] and mask le annulus[1])
; data = fltarr(sc[3],nump)

; writefits, 'masktest.fits', good

 indices=good
 data = fltarr(n_elements(good), sc[3])

                                ; Go through the images, layer by
                                ; layer, and extract the data points
                                ; required here
 for ii=0, sc[3]-1 do begin
    img= cube[*,*,ii]
    data[*,ii] = img[good]
 endfor

return, {data:data, good:good, num:n_elements(good), indices:indices}
end

