function adiklip, cube,  k_klip, trial=trial, target=target, checkpoint = checkpoint, anglemask=anglemask, $
distmask=distmask, posang=posang, diam=diam,pixelscale=pixelscale, wl=wl,angsep=angsep,hyper=hyper,$
obj=obj,nrings=nrings,wr =wr,n_ang =n_ang, num=num, anglemax=anglemax,annmode_inout=annmode_inout,spot_radius=spot_radius,rho=rho,phi=phi

newline = string(10B)

;Function to create PSF subtracted image using KLIP algorithm
;The target image is currently the first image in the PSF cube
;STRUCTURE - Output from search_zones.pro
;IM_SIZE - Dimension of square subarray
;K_KLIP - Number of KL Basis Vectors to retain

; diam - telescope diameter in m 
; wl - wavelength of the observations
; pixelscale - in 1/arcsec


; hyper - if set, KLIP will only process four round patches centered
;         on the planets, but not the rest of the images
; =======================

 
if keyword_set(hyper) then begin
         print, '--[ ADIKLIP    ]-- STARTING HYPERKLIP'
         if keyword_set(trial) then print, 'Trial:', trial
nimages = (size(cube))[3]
;nspots = 1  ; Number of spots  = number of planets
;wr = 7	;MagAO
;wr = 12    ;SPHERE
			; radius of the spot in pixels
;phi=[17.]*!DTOR & rho=[0.7]/pixelscale
;if num eq 131399 then begin nspots=8 & phi=[196.5,16.5,106.5,286.5,196.5+45.,16.5+45.,106.5+45.,286.5+25.,196.5+25.,16.5+25.,106.5+25.,286.5+25.]*!DTOR & rho=[0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85]/pixelscale & endif
;num=1;31399
;if num eq 131399 then begin

	;nspots=6. ;MagAO
	;nspots=7. ;SPHERE
	;nspots=1.
	
	nspots=n_elements(rho)
	;phi=fltarr(nspots)
	;rho=fltarr(nspots*4.)
	;rho[*]=0.79/pixelscale
	;for ii=0.,nspots-1. do phi[ii]=( 196.7+(360.* ii / nspots) )*!DTOR
	;phi=[phi,phi,phi,phi]	;magao
	
	;for ii=0.,nspots-1. do phi[ii]=( 193.+(360.* ii / nspots) )*!DTOR
	
	;rho=rho/pixelscale ;rho is input in arcsec
	;phi=phi*!DTOR
	
	wr=spot_radius
	
	;phi=phi-90.
	
	;PHI=194.*!DTOR
	
;	phi=120.*!DTOR
	;;phi=194.+45.
	;phi=194.
	;phi=194.+20.
	;phi=194.-36.
	;phi=[phi,phi+90.,phi-90.,phi-180.]*!DTOR	;sphere
	
	;rho[0:nspots-1]=0.822/pixelscale ;MagAO
	;rho[0:nspots-1]=0.835/pixelscale	;SPHERE

	
	;rho[nspots:2*nspots-1]=0.55/pixelscale
	;rho[2*nspots:3*nspots-1]=0.4/pixelscale
	;rho[3*nspots:4*nspots-1]=1./pixelscale

	;nspots=1;nspots*4


;endif




final_arr = reform(cube[*,*,0]) ; final array same size as sub array
replicate_inplace, final_arr, !values.f_nan


;Build final image segment by segment:
for ii=0, nspots-1 do begin
 print, '--[ ADIKLIP  ]-- Processing Spot #',ii
      ; What pos angle difference translates to 1.5 FWHM in this annulus/ring?

   ;if nspots gt 1 then arcdist = (asin((1.22 * (wl *(10.^(-6.)))/diam/!DTOR *3600./pixelscale ) / (rho[ii])) )/!DTOR ; Gives the minimum posang difference
	;if nspots le 1 then arcdist = (asin((1.22 * (wl *(10.^(-6.)))/diam/!DTOR *3600./pixelscale ) / (rho)) )/!DTOR ; Gives the minimum posang difference

    if nspots gt 1 then  arcdist= 360.*((wl*(1.0E-6)/diam)*(206265.)/pixelscale ) / (2.*!PI*rho[ii])   $
    	else  arcdist= 360.*((wl*(1.0E-6)/diam)*(206265.)/pixelscale ) / (2.*!PI*rho)                  ; width wr

                                ; required to have a 1.5 FWHM arc
                                ; difference in H-band with a
                                ; diffraction limited 8m telescope at
                                ; the mid-point of the iith ring with
                                ; width wr
                                print, rho;[ii]
    print, 'Arcdist constraint = ', ANGSEP*arcdist   
    arcdist=reform(arcdist)                         
 
                                ; AD March 3, 2015

;first, we get the data from our target segment (the one we will
;project onto KL basis

  ; Calculate the planets' position angle consider field rotation

 zones = getpatch(cube, rho[ii], phi[ii]- 1.*mean(posang)*!DTOR,wr)  ; Get patch corresponding to the planet ii

;        zones = getzone(cube, annulus, angles)  ; Returns a nz x npix 2D array with the pixels in the given zone across the cube

      target_seg = reform(zones.data[*,target]) 

      frames = findgen(nimages)
    ;  refs = indgen(nimages)
     ; print, size(arcdist)
     ; print, refs
     ; print, abs(posang-posang[target])
     ; print, arcdist*angsep
      thresh=arcdist*angsep
     ; print, where(abs(posang-posang[target]) gt thresh[0])
;refs=indgen(nimages)
	;  nottarget  =  where( abs(posang[frames]-posang[target])  ge ANGSEP*abs(arcdist)   ) ; AD
	;print, 'Building reference with frames: ', abs(posang[frames]-posang[target])
    ;  nottarget  =  refs[where( abs(posang[frames]-posang[target])  gt ANGSEP*(arcdist) )]; and abs(posang[frames]-posang[target])  LE anglemax ) ; AD
  	nottarget=where(abs(posang-posang[target]) gt thresh[0])
  	print, 'Not target = ', n_elements(nottarget)
  
  
      if N_elements(nottarget) lt 4 then begin
           nottarget =  where(findgen(nimages) ne target)
          print, ' -[ ADIKLIP ]- Less than 4 frames avialable for PCA - Removing angular distance constraint'
        endif
     
      data_seg = transpose(zones.data[*,nottarget])
 
      dense_target = transpose(target_seg - mean(target_seg))
;      dense_target = target_seg[where(finite(target_seg) eq 1,/null)] ; Remove NaN   

;      dense_target = dense_target - total(dense_target)/n_elements(dense_target)  ; subtract mean
;      dense_target = transpose(dense_target)  ; 1XN

;      size_target = n_elements(dense_target)
;      size_mask = n_elements(zones.indices)

;      dense_data = reduce_matrix(data_seg)     ; Condense Matrix (rid NaN)

      dense_data = data_seg
      size_data = (size(dense_data))[2]

      kl_basis = get_klip_basis(dense_data, k_klip) ; CREATE BASIS!

;last, we will project our target image onto our KL basis

      synthetic_1d = dense_target # matrix_multiply(kl_basis,kl_basis,/atranspose)  ; create psf reconstruction
      
      synthetic_1d = dense_target - synthetic_1d ; subtract psf reconstruction     
                 
      final_arr[zones.indices] = synthetic_1d ; Build Final Image  -- AD: Replaced loop with single command


   
endfor

   


endif 

if keyword_set(annmode_inout) then begin


print, '--[ ADIKLIP    ]-- STARTING ANNULAR KLIP'
if keyword_set(trial) then print, 'Trial:', trial
;seperate structure elements
;ringseg = structure.ringseg
;mask_index = structure.mask_index
;structure = 0.0  ; free memory


;store number of images, rings, and angular segments in ringseg
;nimages = (size(ringseg))[1]
;nrings = (size(ringseg))[2]
;n_ang = (size(ringseg))[3]
;data = (size(ringseg))(4)

nimages = (size(cube))[3]

;Set dimensions of single frame
n_cols = (size(cube))[1]
n_rows = n_cols

final_arr = replicate(!values.f_nan, n_cols, n_rows); final array same size as sub array

;full_lib = fltarr(nimages-1,nrings,n_ang,data)  ; Define array for psf lib

;for aa = 0, nimages-2 do begin
;   full_lib[aa,*,*,*] = ringseg[aa+1,*,*,*]  ; Remove target image from psf lib
;endfor

;  full_lib[0:nimages-2,*,*,*] = ringseg[1:nimages-1,*,*,*]  ; Replaced the above loop with this single line

;store number of images in psf lib
;nimages = (size(full_lib))[1]

progress_counter = 0.0          ; initialize progress counter
total_elements = nrings*n_ang   ; also used for progress update

;Build final image segment by segment:
;for ii=0, nrings-1 do begin

      ; What pos angle difference translates to 1.5 FWHM in this annulus/ring?

   arcdist = (asin((1.22 * (wl * (10.^(-6.))  )/diam/!DTOR *3600./pixelscale ) / (annmode_inout[0]+(abs(annmode_inout[1]-annmode_inout[0] )/2.)) ))/!DTOR ; Gives the minimum posang difference
                                ; required to have a 1.5 FWHM arc
                                ; difference in H-band with a
                                ; diffraction limited 8m telescope at
                                ; the mid-point of the iith ring with
   arcdist= 360.*((wl*(1.0E-6)/diam)*(206265.)/pixelscale ) / (2.*!PI*(annmode_inout[0]+(abs(annmode_inout[1]-annmode_inout[0] )/2.)))                          ; width wr
 	print, 'Arcsist constraint:',arcdist*angsep
                         
 
                                ; AD March 3, 2015
;  print, 'Ring #',ii,' with min angular distance of :',arcdist
   for jj=0, n_ang-1 do begin


;first, we get the data from our target segment (the one we will
;project onto KL basis

       annulus = [ annmode_inout[0],annmode_inout[1] ]
     ;  print, annulus
       angles   = [ jj*(360./n_ang) , (jj+1)*(360./n_ang) ]
       if max(angles) gt 360 then stop

              ;  print, annulus
              ;  print, angles

        zones = getzone(cube, annulus, angles)  ; Returns a nz x npix 2D array with the pixels in the given zone across the cube

      target_seg = reform(zones.data[*,target]) 

;      nottarget =  where(findgen(nimages) ne target) ; AD
                                ; Changing this to exclude frames with
                                ; relative pos angles < 1.5 FWHM at
                                ; the given distance
      frames = findgen(nimages)

      nottarget  =  where( abs(posang[frames]-posang[target])  ge ANGSEP*abs(arcdist) and abs(posang[frames]-posang[target])  LE anglemax  ) ; AD

      ;nottarget  =  where( abs(posang(frames)-posang(target))  ge ANGSEP*abs(arcdist) ) ; AD
          	if jj eq 0 then print, 'Frames meeting criteria:', n_elements(nottarget)
	if jj eq 0 then print, 'Annmode radius = ',(annmode_inout[0]+(abs(annmode_inout[1]-annmode_inout[0] )/2.))
      if N_elements(nottarget) lt 4 then begin
           nottarget =  where(findgen(nimages) ne target)
           
          print, ' -[ ADIKLIP ]- Less than 4 frames avialable for PCA on annulus –  - Removing angular distance constraint'

        endif

      data_seg = transpose(zones.data[*,nottarget])
 
      dense_target = target_seg[where(finite(target_seg) eq 1,/null)] ; Remove NaN   

      if dense_target eq !NULL then begin
         
         print, 'ERROR: EMPTY SEARCH ZONE'
         print,'MAKE SURE ALL SEARCH ZONES CONTAIN >1 PIXELS'

      endif

      dense_target -= total(dense_target)/n_elements(dense_target)  ; subtract mean
      dense_target = transpose(dense_target)  ; 1XN
;      index = where(finite(mask_index[ii,jj,*]) eq 1,/null)  ; get segment index

      size_target = n_elements(dense_target)
      size_mask = n_elements(zones.indices)
      if size_target ne size_mask then print, 'MASK NOT SAME SIZE AS DATA'


;next, we will create our KL basis

     
      ;store the segments of interest
      ;data_seg = full_lib[*,ii,jj,*]   ; create psf lib for segment

      dense_data = reduce_matrix(data_seg)     ; Condense Matrix (rid NaN)

      size_data = (size(dense_data))[2]
      
      ;If the segment contains only a single pixel then dense_data is 1d matrix
      if size_target eq 1 then begin

         print, '(WARNING!) Some Search Zones Contain Only 1 Pixel...'
         
         size_data = 1

      endif

      if size_target ne size_data then print, 'TARGET SEG NOT SAME SIZE AS DATA SEG'

      kl_basis = get_klip_basis(dense_data, k_klip) ; CREATE BASIS!

;last, we will project our target image onto our KL basis


      synthetic_1d = dense_target # matrix_multiply(kl_basis,kl_basis,/atranspose)  ; create psf reconstruction
      
      synthetic_1d = dense_target - synthetic_1d ; subtract psf reconstruction     
      
      final_arr[zones.indices] = synthetic_1d ; Build Final Image  -- AD: Replaced loop with single command

      progress_counter += 1 ; update progress

   endfor
   
;endfor

endif 

if not keyword_set(hyper) and not keyword_set(annmode_inout) then begin

print, '--[ ADIKLIP    ]-- STARTING REGULAR KLIP'

if keyword_set(trial) then print, 'Trial:', trial
;seperate structure elements
;ringseg = structure.ringseg
;mask_index = structure.mask_index
;structure = 0.0  ; free memory


;store number of images, rings, and angular segments in ringseg
;nimages = (size(ringseg))[1]
;nrings = (size(ringseg))[2]
;n_ang = (size(ringseg))[3]
;data = (size(ringseg))(4)

nimages = (size(cube))[3]

;Set dimensions of single frame
n_cols = (size(cube))[1]
n_rows = n_cols

final_arr = replicate(!values.f_nan, n_cols, n_rows); final array same size as sub array

;full_lib = fltarr(nimages-1,nrings,n_ang,data)  ; Define array for psf lib

;for aa = 0, nimages-2 do begin
;   full_lib[aa,*,*,*] = ringseg[aa+1,*,*,*]  ; Remove target image from psf lib
;endfor

;  full_lib[0:nimages-2,*,*,*] = ringseg[1:nimages-1,*,*,*]  ; Replaced the above loop with this single line

;store number of images in psf lib
;nimages = (size(full_lib))[1]

progress_counter = 0.0          ; initialize progress counter
total_elements = nrings*n_ang   ; also used for progress update

;Build final image segment by segment:
for ii=0, nrings-1 do begin

      ; What pos angle difference translates to 1.5 FWHM in this annulus/ring?

   arcdist = (asin((1.22 * (wl * (10.^(-6.)))/diam/!DTOR *3600./pixelscale ) / ((float(ii)+0.5)*wr )) )/!DTOR ; Gives the minimum posang difference
                                ; required to have a 1.5 FWHM arc
                                ; difference in H-band with a
                                ; diffraction limited 8m telescope at
                                ; the mid-point of the iith ring with
                                ; width wr
                                
 
                                ; AD March 3, 2015
;  print, 'Ring #',ii,' with min angular distance of :',arcdist
   annulus = [ ii*wr,(ii+1)*wr ]
   frames = findgen(nimages)
   nottarget  =  where( abs(posang[frames]-posang[target])  ge ANGSEP*abs(arcdist) and abs(posang[frames]-posang[target])  LE anglemax  ) ; AD
   for jj=0, n_ang-1 do begin


;first, we get the data from our target segment (the one we will
;project onto KL basis

       angles   = [ jj*(360./n_ang) , (jj+1)*(360./n_ang) ]
       if max(angles) gt 360 then stop

              ;  print, annulus
              ;  print, angles

      zones = getzone(cube, annulus, angles)  ; Returns a nz x npix 2D array with the pixels in the given zone across the cube

      target_seg = reform(zones.data[*,target]) 

;      nottarget =  where(findgen(nimages) ne target) ; AD
                                ; Changing this to exclude frames with
                                ; relative pos angles < 1.5 FWHM at
                                ; the given distance

      ;nottarget  =  where( abs(posang(frames)-posang(target))  ge ANGSEP*abs(arcdist) ) ; AD
  
      if N_elements(nottarget) lt 4 then begin
           nottarget =  where(findgen(nimages) ne target)
          if jj eq 0 then print, ' -[ ADIKLIP ]- Less than 4 frames avialable for PCA on annulus – ',ii,' - Removing angular distance constraint'
        if jj eq 0 then print, 'Arcsist constraint:',arcdist
        	;print, posang
      endif

      data_seg = transpose(zones.data[*,nottarget])
 
      dense_target = target_seg[where(finite(target_seg) eq 1,/null)] ; Remove NaN   

      if dense_target eq !NULL then begin
         
         print, 'ERROR: EMPTY SEARCH ZONE'
         print,'MAKE SURE ALL SEARCH ZONES CONTAIN >1 PIXELS'

      endif

      dense_target -= total(dense_target)/n_elements(dense_target)  ; subtract mean
      dense_target = transpose(dense_target)  ; 1XN
;      index = where(finite(mask_index[ii,jj,*]) eq 1,/null)  ; get segment index

      size_target = n_elements(dense_target)
      size_mask = n_elements(zones.indices)
      if size_target ne size_mask then print, 'MASK NOT SAME SIZE AS DATA'


;next, we will create our KL basis

     
      ;store the segments of interest
      ;data_seg = full_lib[*,ii,jj,*]   ; create psf lib for segment

      dense_data = reduce_matrix(data_seg)     ; Condense Matrix (rid NaN)

      size_data = (size(dense_data))[2]
      
      ;If the segment contains only a single pixel then dense_data is 1d matrix
      if size_target eq 1 then begin

         print, '(WARNING!) Some Search Zones Contain Only 1 Pixel...'
         
         size_data = 1

      endif

      if size_target ne size_data then print, 'TARGET SEG NOT SAME SIZE AS DATA SEG'

      kl_basis = get_klip_basis(dense_data, k_klip) ; CREATE BASIS!

;last, we will project our target image onto our KL basis


      synthetic_1d = dense_target # matrix_multiply(kl_basis,kl_basis,/atranspose)  ; create psf reconstruction
      
      synthetic_1d = dense_target - synthetic_1d ; subtract psf reconstruction     
      
      final_arr[zones.indices] = synthetic_1d ; Build Final Image  -- AD: Replaced loop with single command

      progress_counter += 1 ; update progress

   endfor
   
endfor

endif
;===================

print, '--[ ADIKLIP    ]-- COMPLETE' + newline

return, final_arr
END
