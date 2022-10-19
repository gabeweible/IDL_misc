; 12.12.99: implementation for rectangular images

; 14.1.99: bug for angle 90 degrees removed

; removes sinusodal pattern on _square_ images
; d - image or cube of images
; angle - inclination of pattern with respect to x axis
; clip_level - clip level above which image is clipped in median calculation
; for MAGIC images, it is recommended to clean the four subframes individually
; case for 0 and 90 degrees do not require rotation and are treated separately

function destripe,d,angle,clip_level=clip_level,comment=comment,no_fit=no_fit,$
        nodisp=nodisp,fraction=fraction

; on_error,2

; median should be computed on at least a fraction of all pixels!
if not keyword_set(fraction) then fraction=0.5

n=size(d)	; get image dimension
nx=n(1)
nx2=nx/2
ny=n(2)
ny2=ny/2
if (n(0) eq 2) then n3=1 else n3=n(3)
h=d		; h is temporary cube, holds cleaned frames

r=angle/!radeg
m=[[cos(r),sin(r)],[-sin(r),cos(r)]]
nxx=ceil(nx*max([(m#[1,1])[0],(m#[1,-1])[0],(m#[-1,-1])[0],(m#[-1,1])[0]]))
nyy=ceil(ny*max([(m#[1,1])[1],(m#[1,-1])[1],(m#[-1,-1])[1],(m#[-1,1])[1]]))

hhh=fltarr(nxx,nyy)	; holds rotated single image
mask=bytarr(nxx,nyy)	; mask for median computation
mask((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1)=1
if keyword_set(clip_level) then $
 if (n3 gt 1) then mask((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1)= $
                        mask((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1) $
                         and (total(d,3)/n3 lt clip_level) $
 else mask((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1)= $
        mask((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1) and (d lt clip_level)
if not keyword_set(nodisp) then tvscl,mask((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1)
xp=findgen(nxx)
yp=findgen(nyy)
case angle of		; check first for obvious cases
0: begin
;	weight=medcol(mask,mask,1)
        weight=total(mask,1)
	for i=0,n3-1 do begin
	 hhh((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1)=d(*,*,i)
   	 x=medcol(hhh,mask,1)		; compute median rowwise
 	 index=where(weight gt fraction*nx,cnt)
; are there enough points for interpolation?
         if cnt lt 3 then begin
                message,'No action possible: lower clip_level/fraction!'
         end
         if keyword_set(no_fit) then st=replicate(1.,nxx)#(x*(weight gt fraction*nx)) else begin
                dummy=SPL_INIT(yp[index], x[index])
                xf=SPL_INTERP(yp[index], x[index], dummy, yp)
	        st=replicate(1.,nxx)#xf
         end
	 h(*,*,i)=d(*,*,i)-st((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1)
	 if keyword_set(comment) then print,'loop ',i,' finished'
	endfor	
   end
90: begin
; 	weight=medcol(mask,mask,2)		
        weight=total(mask,2)
	for i=0,n3-1 do begin
		hhh[where(hhh eq 0)]=!values.f_nan
		hhh[where(hhh eq 0.)]=!values.f_nan

	 hhh((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1)=d(*,*,i)
   	 x=medcol(hhh,mask,2)		; compute median columnwise
         index=where(weight gt fraction*ny,cnt)
; are there enough points for interpolation?
         if cnt lt 3 then begin
                message,'No action possible: lower clip_level/fraction!'
         end
	 if keyword_set(no_fit) then st=x*(weight gt fraction*ny)#replicate(1.,nyy) else begin
                xf=INTERPOL(x[index], xp[index], xp)
                st=xf#replicate(1.,nyy)				 
         end
	 h(*,*,i)=d(*,*,i)-st((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1)
	 if keyword_set(comment) then print,'loop ',i,' finished'
	endfor
   end
else: begin
 mask=rot(mask,angle-90,/cubic)	; rotate to vertical orientation		; y vector
; weight=medcol(mask,mask,2)
 weight=total(mask,2)
 for i=0,n3-1 do begin		; treat single images
  hhh((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1)=d(*,*,i)
  hhh=rot(hhh,angle-90,cubic=-0.5)	; rotate to vertical orientation
  x=medcol(hhh,mask,2)		; compute median columnwise
  index=where(weight gt fraction*ny,cnt)
; are there enough points for interpolation?
  if cnt lt 3 then begin
        message,'No action possible: lower clip_level/fraction!'
  end
  if keyword_set(no_fit) then st=x*(weight gt fraction*ny)#replicate(1.,nyy) else begin
  	  dummy=poly_fit(yp(index),x(index),4,xx)	; this is experimental
  	  xf=x					; try to avoid introducing spurious features
  	  xf(index)=xx				; on spatial scales larger than ripples
          st=(x-xf)#replicate(1.,nyy)				 
  end
  st=rot(st,90-angle,cubic=-0.5)	; create image and rotate back
  h(*,*,i)=d(*,*,i)-st((nxx-nx)/2:(nxx-nx)/2+nx-1,(nyy-ny)/2:(nyy-ny)/2+ny-1)
  if keyword_set(comment) then print,'loop ',i,' finished'
 endfor
end
endcase
return,h
end


