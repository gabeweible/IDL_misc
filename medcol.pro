; computes median on rows or columns of a _square_ image (x_or_y=1,2)
; mask can be used to reject image areas from computation
; IDL median seems not to give an good estimate even when mask is applied
; therefore, we should perhaps switch to mmm

function medcol,d,mask,x_or_y


on_error,2
n=size(d)
if (x_or_y eq 1) then begin
 hhh=fltarr(n(2))
 for i=0,n(2)-1 do begin
  index=where(mask(*,i) gt 0,count)
  if (count eq 0) then hhh(i)=0 else hhh(i)=median(d(index,i),/even)
;   if (count eq 0) then hhh(i)=0 else begin
;    mmm,d(index,i),m,s
;    hhh(i)=m
;   end
 end
end
if (x_or_y eq 2) then begin
 hhh=fltarr(n(1))
 for i=0,n(1)-1 do begin
  index=where(mask(i,*) gt 0,count)
  if (count eq 0) then hhh(i)=0 else hhh(i)=median(d(i,index),/even)
;   if (count eq 0) then hhh(i)=0 else begin
;    mmm,d(i,index),m,s
;    hhh(i)=m
;   end

 end
end
return,hhh
end


