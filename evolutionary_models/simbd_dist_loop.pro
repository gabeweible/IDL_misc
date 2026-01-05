function simbd_dist_loop, mjup_arr=mjup_arr, age=age,$
 sigage=sigage, d=d, sigd=sigd, filter=filter,$
 min_eff=min_eff, n_trials=n_trials
 
 compile_opt IDL2
 newline = string(10B)
 
 mags = 0
 print, (size(mjup_arr))[1], ' masses'
 
 for ii=0, (size(mjup_arr))[1]-1 do begin
 
 	mass_ii = mjup_arr[ii]
 	
 	    print, 'Mass index ', ii, ' with mass ', mass_ii, ' Mjup'
 	
 	mag_ii = simbd_dist(mass=mass_ii, age=age, sigage=sigage, d=d, sigd=sigd, filter=filter,$
 	 min_eff=min_eff, n_trials=n_trials)
 	 
 	 mag_ii = mean(mag_ii)
 	 
 	 print, filter, ' mean mag: ', mag_ii
 	 
 	if ii eq 0 then begin
 		mags = mag_ii
 	endif else begin
 		mags = [mags, mag_ii]
 	endelse
 	
 	print, 'mags: ', mags, newline
 	
 endfor
 
 return, mags
 
 end