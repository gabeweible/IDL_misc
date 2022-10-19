;Takes a cube with different image rotation angles and derotates and combines,
;performing a noise-weighted combination of the pixel values for the final image

;making this routine a drop-in for median or mean combination...

;assuming input has already been derotated

;possible to add running variance

; Written by Kevin Wagner
; Based off of Bottom et al. 2017: https://arxiv.org/pdf/1711.09119.pdf

function nw_ang_comb, in_cube, angles

;copy cube so as to leave original input alone
adi_cube=in_cube

;first go back to pupil stabilization
for ii=0, (size(adi_cube))(3)-1 do adi_cube[*,*,ii]=rot(adi_cube[*,*,ii],angles[ii],/interp)

;need to replace NaN for the rest to work correctly
adi_cube[where(finite(adi_cube) lt 1)]=0.

;compute variance of each pixel with time
var=variance(adi_cube,dim=3,/double)


;copy into 3d cube
var_cube=adi_cube
for nn=0,(size(adi_cube))(3)-1 do var_cube[*,*,nn]=var

;derotate cubes
for ii=0, (size(adi_cube))(3)-1 do begin
	adi_cube[*,*,ii]=rot(adi_cube[*,*,ii],-angles[ii],/interp)
	var_cube[*,*,ii]=rot(var_cube[*,*,ii],-angles[ii],/interp)
endfor


;generate an array of ones for division purposes
one_cube=var_cube
one_cube[*]=1.0

;setting up images, a la step 6 in Bottom et al. 2017
im1=one_cube/((size(adi_cube))(3)*mean(one_cube/var_cube,dim=3))
im2=(size(adi_cube))(3)*mean(adi_cube/var_cube,dim=3)


;perform element wise multiplication
nwadi=im1*im2

;return the combined image
return, nwadi

end
