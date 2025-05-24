function bin_frames, frame_folder, angles, szbin, bin_type=bin_type

compile_opt IDL2
newline=string(10B)

if not keyword_set(bin_type) then bin_type='mean'

return_list = list()

if bin_type eq 'mean' then begin
	print, 'Searching for FITS files in', frame_folder, '...'
	files = FILE_SEARCH(frame_folder, '*.fits', COUNT=filecount)
	print, 'Found ', filecount, ' FITS files!'

	; Use lists to take advantage of faster appending over arrays (we'll convert back to arrays later)
	print, 'Initializing...'

	; Initialize an empty coadd frame and our angle for the coadd as 0
	coadd_frame = readfits(files[0])
	replicate_inplace, coadd_frame, 0.
	coadd_angle = 0.
	
	return_cube = list()
	return_angles = list()

	; k counts what image we're on, which is used to know when to add the coadded frame to the cube
	k = 1
	; Start looping through each image
	print, 'Initialized successfully, beginning cube creation loop...'
	for ii = 0, filecount-1 do begin
		print, 'File index', ii, '/', filecount-1
		print, 'Reading in:', files[ii], newline
		frame = readfits(files[ii])
		angle = angles[ii]
		
		coadd_frame += frame; add to the big coadd
		coadd_angle += angle
		
		; coadd group complete
		if k mod szbin eq 0 then begin
			; bin angle group and reset
			coadd_angle *= 1. / szbin
			return_angles.Add, coadd_angle, /no_copy
			coadd_angle = 0.
			
			;Add the mean frame to the return cube and reset the frame
			coadd_frame *= 1. / szbin
			return_cube.Add, [[coadd_frame]], /no_copy
			replicate_inplace, coadd_frame, 0.
		endif; coadd group if
		
		print, 'Number of Frames in Cube: ', n_elements(return_cube)
		k += 1; increment counter
	   
	endfor; loop over frames for
	
	; convert our lists to arrays
	print, 'Converting to arrays...'
	return_cube = return_cube.toArray(/TRANSPOSE, /NO_COPY)
	return_angles = return_angles.toArray(/NO_COPY)
	
	return_list.Add, return_cube, /no_copy
	return_list.Add, return_angles, /no_copy

endif; 'mean' binning if

return, return_list

end; end the program