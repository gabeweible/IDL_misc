;\tiny
;\begin{verbatim}

function angmask, img, xp, yp, degree=degree,dmin=dmin,dmax=dmax

; Creates Angular Mask:
; The value of each point describes its angle from the center


; If dmin and/or dmax is defined, it returns a mask, which
; has values of 1. between angles in [dmin,dmax] range else
; 0.
; AD 2003
;---------------------------





msk=img
sx=(size(msk))(1)
sy=(size(msk))(2)

if not keyword_set(xp) then xp = ((sx-1)/2.)
if not keyword_set(yp) then yp = ((sy-1)/2.)

for xx=0, sx-1 do begin
	for yy=0, sy-1 do begin
		;msk(xx,yy)= atan( (0.9999999*xx-xp)/(0.99999*yy-yp) )
		msk(xx,yy)=atan(yy-yp, xx-xp)
		;if xx lt xp and yy gt yp then msk(xx,yy)=msk(xx,yy)+2*!Pi
		;if yy lt yp then msk(xx,yy)=msk(xx,yy)+!Pi
		
		
	endfor
endfor
msk=msk+!Pi

if keyword_set(degree) then msk = 180. * msk /!Pi

if keyword_set(dmax) or keyword_set(dmin) then  begin
	msk = 180. * msk /!Pi
	;print, 'Returning Angle Mask'
	if not keyword_set(dmax) then dmax=360.
	if not keyword_set(dmin) then dmin=0.
	
	dummy=msk
	dummy(*,*)=0.
	dummy(where(msk gt dmin and msk lt dmax))=1. 
	msk=dummy
endif	


return, msk

end

;\end{verbatim}
;\normalsize
