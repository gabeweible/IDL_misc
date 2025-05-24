;+
; NAME:
;  medianw
; PURPOSE: (one line)
;  minimize the weighted average deviation.
; DESCRIPTION:
;  This routine estimates the average data value, xmed, such that
;  the average deviation, total(weight * abs(x-xmed) ), is minimized.
;  For equally-weighted points, this is the median.
;
;  The statistics are robust in that the result is insensitive
;  to outliers.  
; CATEGORY:
;  Statistics
; CALLING SEQUENCE:
;  mean = medianw(x, w)
; INPUTS:
;  x   - Input data to be analyzed.
; OPTIONAL INPUT PARAMETERS:
;  w - weights.  Assumed equally weighted if not passed.
; INPUT KEYWORD PARAMETERS:
; OUTPUT KEYWORD PARAMETERS:
; OUTPUTS:
; COMMON BLOCKS:
;  None.
; SIDE EFFECTS:
;  None.
; RESTRICTIONS:
;  None.
; PROCEDURE:
;
;
; Minimizing f(xmed) = total(weight * abs(x-xmed)) is equivalent to finding
; df(xmed)/d xmed = 0 = -total(weight * sgn(x-xmed)).  
; If x is sorted, and xmed=x[imed], this can be written as
; 0 = total(weight[0:imed-1] * sgn(x[0:imed-1]-xmed)) 
;   + total(weight[imed+1:n-1] * sgn(x[imed+1:n-1]-xmed))
; or
; 0 = - total(weight[0:imed-1]) + total(weight[imed+1:n-1])
;
; or
; total(weight[0:imed-1]) = total(weight[imed+1:n-1])
;
;
; MODIFICATION HISTORY:
;  Written by Leslie A. Young, Soutwest Research Institute, 2002 Jul 23.
;  Modified LAY Oct 2002.  Better medianwTEST; 
;            Avoid eq for comparing floats in test for one median value
;            (equivalent to an odd number of equally weighted values)
;            Average two neighbors if xmed is not one of the listed x's
;            (equivalent to an even number of equally weighted values)
; Modified LAY Oct 2002.  Changed name to medianw; 
;-
function medianw, x, w

compile_opt IDL2

  n = n_elements(x)
  noweight = (n_params() lt 2)
  if noweight then  w = replicate(1.,n)

  if noweight then begin
    med = median(x, /EVEN, /double)
  end else begin
    s = sort(x)
    xx = float(x[s])
    ww = w[s]/total(w, /double, /nan)
    i = 0
    totw = 0.
    while totw lt (1-(totw+ww[i])) and i lt n-1 do begin
      totw = totw + ww[i]
      i = i + 1
    end
    if (abs( totw - (1-(totw+ww[i])) ) lt ww[i]*1e-6) then begin
      med = xx[i]
    end else begin
      med = ( 0.5*xx[i] + 0.5*xx[i-1] )
    end
  end
  
  return, med
  
end

pro medianwTEST
  x = [1,2,3,4,5]
  marr = findgen(31)*0.1+1.5
  ad = fltarr(31)
  test = 0
  print, 'medianwTEST test #', test
  w = [1,1,1,1,1]
  for i=0,30 do ad[i]=total(abs(x-marr[i])*w, /double, /nan)
  print, 'x', x, format='(A14,5(F6.1))'
  print, 'w', w, format='(A14,5(F6.1))'
  print, 'medianw', medianw(x, w), format='(A14,F6.1)'
  print, 'sample xmed', marr, format='(A14,31(F6.1))'
  print, 'abs deviation', ad, format='(A14,31(F6.1))'
  test = test+1
  
  print, 'medianwTEST test #', test
  w = [1,2,5,1,2]
  for i=0,30 do ad[i]=total(abs(x-marr[i])*w)
  print, 'x', x, format='(A14,5(F6.1))'
  print, 'w', w, format='(A14,5(F6.1))'
  print, 'medianw', medianw(x, w), format='(A14,F6.1)'
  print, 'sample xmed', marr, format='(A14,31(F6.1))'
  print, 'abs deviation', ad, format='(A14,31(F6.1))'
  test = test+1
  
  print, 'medianwTEST test #', test
  w = [2,1,1,1,2]
  for i=0,30 do ad[i]=total(abs(x-marr[i])*w)
  print, 'x', x, format='(A14,5(F6.1))'
  print, 'w', w, format='(A14,5(F6.1))'
  print, 'medianw', medianw(x, w), format='(A14,F6.1)'
  print, 'sample xmed', marr, format='(A14,31(F6.1))'
  print, 'abs deviation', ad, format='(A14,31(F6.1))'
  test = test+1
  
  
  print, 'medianwTEST test #', test
  w = [2,2,1,1,2]
  for i=0,30 do ad[i]=total(abs(x-marr[i])*w)
  print, 'x', x, format='(A14,5(F6.1))'
  print, 'w', w, format='(A14,5(F6.1))'
  print, 'medianw', medianw(x, w), format='(A14,F6.1)'
  print, 'sample xmed', marr, format='(A14,31(F6.1))'
  print, 'abs deviation', ad, format='(A14,31(F6.1))'
  test = test+1
  
  print, 'medianwTEST test #', test
  w = -[2,2,1,1,2]
  for i=0,30 do ad[i]=total(abs(x-marr[i])*w)
  print, 'x', x, format='(A14,5(F6.1))'
  print, 'w', w, format='(A14,5(F6.1))'
  print, 'medianw', medianw(x, w), format='(A14,F6.1)'
  print, 'sample xmed', marr, format='(A14,31(F6.1))'
  print, 'abs deviation', ad, format='(A14,31(F6.1))'
  test = test+1
  
  
end
