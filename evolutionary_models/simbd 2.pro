function simbd, mass=mass,age=age,filter=filter, debug=debug, error=error, stp=stp, usepreinterp=usepreinterp,reload=reload,mdusty=mdusty
  
;mass in solar masses.


; Edited by KW on August 11 2018

; Confirmed that this reproduces the right values for non-interpolated
; ages - AD June 26, 2008

;====================================================
; Script outputs the magnitude in the given filter 
; of the simulated BD/EGP atmosphere
; with the given age, mass, etc.
; Uses data from Baraffe et al. 2003
;
; AD2004

; Rewritten to allow interpolation for arbitrary mass and age (within the grid)
; D. Apai, March 7, 2007, Steward Observatory 
;====================================================
;
; Needs data files in the ~/simbd.data/ directory 

;
;====================================================


;mass=mass*0.0009543	;Mj to Msun, commented out because it apepars in simbd_dist

pwd='/Users/gweible/Library/CloudStorage/OneDrive-UniversityofArizona/research/HII1348/big_data/SimBD/'

if keyword_set(mdusty) then datadir='simbd.data.dusty/' else datadir='simbd.data/'

if keyword_set(mdusty) then columns=14 else columns=13

if not keyword_set(reload) then begin
aiuread_ascii, pwd+datadir+'1Myr.dat', one, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'5Myr.dat', five, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'10Myr.dat', ten, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'50Myr.dat', fifty, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'100Myr.dat', hundred, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'120Myr.dat', hundredtwenty, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'500Myr.dat', fivehundred, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'1Gyr.dat', onegig, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'5Gyr.dat', fivegig, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'10Gyr.dat', tengig, ncol=columns,/quiet


if keyword_set(usepreinterp) then begin
; The interpolated values:
aiuread_ascii, pwd+datadir+'12Myr.dat', twelve, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'20Myr.dat', twenty, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'30Myr.dat', thirty, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'70Myr.dat', seventy, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'200Myr.dat', twohundred, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'680Myr.dat', sixeighty, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'2000Myr.dat', twogig, ncol=columns,/quiet
aiuread_ascii, pwd+datadir+'3000Myr.dat', threegig,ncol=columns, /quiet
aiuread_ascii, pwd+datadir+'4000Myr.dat', fourgig,ncol=columns, /quiet

ndata = 19

endif else ndata=10.

singleage = { Age:age, data:one}

simdata = replicate(singleage,ndata) 
simdata(0).data = one
simdata(0).age  = 1.
simdata(1).data = five
simdata(1).age  = 5.
simdata(2).data = ten
simdata(2).age  = 10.
simdata(3).data = fifty
simdata(3).age  = 50.
simdata(4).data = hundred
simdata(4).age  = 100.
simdata(5).data = hundredtwenty  
simdata(5).age  = 120.
simdata(6).data = fivehundred
simdata(6).age  = 500.
simdata(7).data = onegig
simdata(7).age  = 1000.
simdata(8).data = fivegig
simdata(8).age  = 5000.
simdata(9).data = tengig
simdata(9).age  = 10000.


if keyword_set(usepreinterp) then begin

; The pre-interpolated:
simdata(10).data = twelve
simdata(10).age  = 12.
simdata(11).data = twenty
simdata(11).age  = 20.

simdata(12).data = thirty
simdata(12).age  = 30.

simdata(13).data = seventy
simdata(13).age  = 70.
simdata(14).data = twohundred
simdata(14).age  = 200.
simdata(15).data = sixeighty
simdata(15).age  = 680.

simdata(16).data = twogig
simdata(16).age  = 2000.

simdata(17).data = threegig
simdata(17).age  = 3000.

simdata(18).data = fourgig
simdata(18).age  = 4000.


endif

if not keyword_set(mdusty) then save, simdata, filename=pwd+'isochrones.dat' else save, simdata, filename=pwd+'dusty-isochrones.dat'

endif else if not keyword_set(mdusty) then restore, pwd+'isochrones.dat' else restore, pwd+'dusty-isochrones.dat'


;filter='Mh'
; ===============================================================================================
    	field=['Mass','Teff','Lum','g','R','Mv','Mr','Mi','Mj','Mh','Mk','Ml','Mm'] 

	if keyword_set(mdusty) then   field=['Mass','Teff','Lum','g','R','li','Mv','Mr','Mi','Mj','Mh','Mk','Ml','Mm'] 

    	sel_filter = where(strmatch(field, filter, /fol) eq 1) ; which column?

; ===============================================================================================
; CREATE DATA CUBE FOR THE GIVEN FILTER
; ===============================================================================================

    ages   = reform(simdata[*].age)
    masses = reform(simdata[0].data[0,*])


;create matrix, remove missing points, and upscale with bilinear interpolation, commented out for reloading purposes
if not keyword_set(reload) then begin
   matrix = reform(simdata[*].data[sel_filter[0],*]) ; Make 2D matrix with masses and age axes
   if keyword_set(mdusty) then matrix[where(matrix eq -1)]=!values.f_nan else matrix[where(matrix eq -1)]=max(matrix)
   matrix=congrid(matrix,n_elements(masses)*50.,n_elements(ages)*50.,/interp)



if not keyword_set(mdusty) then if filter eq 'Mk' then writefits,pwd+'cond-matrix-congrid-K-50.fits',matrix 
if not keyword_set(mdusty) then if filter eq 'Ml' then  writefits,pwd+'cond-matrix-congrid-L-50.fits',matrix
if not keyword_set(mdusty) then if filter eq 'Mh' then  writefits,pwd+'cond-matrix-congrid-H-50.fits',matrix

if keyword_set(mdusty) then if filter eq 'Mk' then writefits,pwd+'dusty-matrix-congrid-K-50-nan.fits',matrix 
if keyword_set(mdusty) then if filter eq 'Ml' then writefits,pwd+'dusty-matrix-congrid-L-50-nan.fits',matrix
if keyword_set(mdusty) then if filter eq 'Mh' then writefits,pwd+'dusty-matrix-congrid-H-50-nan.fits',matrix


endif else begin ;reloading makes things go much faster

if not keyword_set(mdusty) then if filter eq 'Mk' then matrix=readfits(pwd+'cond-matrix-congrid-K-50.fits',/silent) else matrix=readfits(pwd+'cond-matrix-congrid-L-50.fits',/silent)

if keyword_set(mdusty) then if filter eq 'Mk' then matrix=readfits(pwd+'dusty-matrix-congrid-K-50-nan.fits',/silent) else matrix=readfits(pwd+'dusty-matrix-congrid-L-50-nan.fits',/silent)

endelse


;upscale masses and ages; be sure that the upscaling factor is the same as the matrix, and also that the interpolation scheme is the same
   masses=congrid(masses,n_elements(masses)*50.,/interp)
      ages=congrid(ages,n_elements(ages)*50.,/interp)


    error  = '' ;create string to contain error messages should they appear later

; ===============================================================================================
; FIND BRACKETING AGES AND MASSES:
; ===============================================================================================

;set age to nearest item on grid, disables interpolation
;age=nearest_element(age,ages)
;mass=nearest_element(mass,masses)


    ; FIND BRACKETING AGES:
    ; ===========================
     agendx = where(ages le age, cnt)
     if cnt eq 0 then lowerage=-1 else lowerage = max(ages[agendx])
     lageindex = where(ages eq lowerage)
    
     agendx = where(ages ge age, cnt)
     if cnt eq 0 then upperage=-1 else upperage = min(ages[agendx])
     uageindex = where(ages eq upperage)
    
    ; FIND BRACKETING MASSES:
    ; ============================
     massndx = where(masses le mass, cnt)
     if cnt eq 0 then lowermass=-1 else lowermass = max(masses[massndx])
     lmassindex = where(masses eq lowermass)


     massndx = where(masses ge mass, cnt)
     if cnt eq 0 then uppermass=-1 else uppermass = min(masses[massndx])
     umassindex = where(masses eq uppermass)

	if n_elements(lageindex) gt 1 then lageindex=lageindex[0]
	if n_elements(uageindex) gt 1 then uageindex=uageindex[0]


	if n_elements(lmassindex) gt 1 then lmassindex=lmassindex[0]
	if n_elements(umassindex) gt 1 then umassindex=umassindex[0]


    if keyword_set(debug) then  print, 'DEBUG: Bracketing Ages:',lowerage, upperage, ' Bracketing Masses:',lowermass, uppermass

; ===============================================================================================
; INTERPOLATE TO THE GIVEN MASS AND AGE
; ===============================================================================================


;print, reform([lageindex, uageindex, lmassindex, umassindex])
    if min([lageindex, uageindex, lmassindex, umassindex]) gt -1 then begin
    
    	    if min( [ matrix[lmassindex,lageindex],matrix[umassindex, uageindex],$
	    matrix[umassindex,lageindex],matrix[lmassindex,uageindex]]) ge 0 then begin
	    ; Bracketed by four values, all of which have a real magnitude value associated to 
	    ;  => interpolation possible
	
	    if keyword_set(debug) then print, 'DEBUG: Found four valid gridpoints around the target value, interpolating. '

		
	    mag1 = matrix[lmassindex,lageindex]
	    mag2 = matrix[umassindex,lageindex]
	    mag3 = matrix[lmassindex,uageindex]
	    mag4 = matrix[umassindex,uageindex]
  

            agedist = upperage-lowerage

            massdist = uppermass-lowermass
            dists = fltarr(4)

            if agedist eq 0. then nagedistlower = 0. else nagedistlower = (lowerage-age)/agedist
            if agedist eq 0. then nagedistupper = 0. else nagedistupper = (upperage-age)/agedist
            if massdist eq 0. then nmassdistlower = 0. else nmassdistlower = (lowermass-mass)/massdist
            if massdist eq 0. then nmassdistupper = 0. else nmassdistupper = (uppermass-mass)/massdist

	    dists[0] = sqrt( (nagedistlower^2.) + (nmassdistlower^2.))
	    dists[1] = sqrt( (nagedistlower^2.) + (nmassdistupper^2.))
	    dists[2] = sqrt( (nagedistupper^2.) + (nmassdistlower^2.))
	    dists[3] = sqrt( (nagedistupper^2.) + (nmassdistupper^2.))



         

            if keyword_set(debug) then  print, 'DEBUG: Mag1 and Distance:', mag1, dists[0]
	    if keyword_set(debug) then  print, 'DEBUG: Mag2 and Distance:', mag2, dists[1]
	    if keyword_set(debug) then  print, 'DEBUG: Mag3 and Distance:', mag3, dists[2]
	    if keyword_set(debug) then  print, 'DEBUG: Mag4 and Distance:', mag4, dists[3]    

           
            w=fltarr(4) & w=double(w)
    	    w[0] = 1.0/dists[0]
	    w[1] = 1.0/dists[1]
	    w[2] = 1.0/dists[2]
	    w[3] = 1.0/dists[3]

            totw = total(w)
            w = w/totw

	


           if min(dists eq [0,0,0,0]) eq 0 then  $
               magnitude = ( (mag1*w[0]) + (mag2*w[1]) + (mag3*w[2]) + (mag4*w[3]) ) else $
                 magnitude = mean([mag1,mag2,mag3,mag4])


    	    if keyword_set(debug) then print, 'DEBUG: Interpolated value:', magnitude	   	 

    	endif else begin

	    if lmassindex eq -1 then error=error+' - Mass beyond the lower end of the grid -'
	    if umassindex eq -1 then error=error+' - Mass beyond the upper end of the grid -'
	    if lageindex eq -1 then error=error+' - Age beyond the lower end of the grid -'
	    if uageindex eq -1 then error=error+' - Age beyond the upper end of the grid -'
	    
	    
	    ; Bracketed by four values, but at least one magnitude equals -1, i.e. lacks
	    ; real data.
	    ; Cannot interpolate, return NaN - this usually happens at very low masses and old ages
	    magnitude = !values.f_nan
    	if keyword_set(debug) then print, 'DEBUG: Inside the grid, but not all four points have real data - cannot interpolate'

	endelse
	
    endif else begin
    	; Outside or on the edge of the grid
    	; return -1
    	if keyword_set(debug) then print, 'DEBUG: Outside the grid points, cannot interpolate'
    	    magnitude = !values.f_nan
	    if lmassindex eq -1 then error=error+' - Mass beyond the lower end of the grid -'
	    if umassindex eq -1 then error=error+' - Mass beyond the upper end of the grid -'
	    if lageindex eq -1 then error=error+' - Age beyond the lower end of the grid -'
	    if uageindex eq -1 then error=error+' - Age beyond the upper end of the grid -'

    
    endelse
    
    
    

; ===============================================================================================
; OUTPUT RESULT
; ===============================================================================================
if keyword_set(stp) then stop

return, magnitude
end
