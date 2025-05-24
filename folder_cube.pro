pro folder_cube, find_path, save_path=save_path, type=type
; Searches for all .FITS files in find_path (string), then combines them into a cube written to
; save_path (string)
; type is 'list' or 'array'
; 
; One can concatenate to lists much faster in IDL, but they are then converted into arrays in order
; to use writefits.pro, which only accepts arrays as input.

compile_opt IDL2
newline = string(10B)

starttime=systime(/JULIAN)
; Search for the raw images in the specified data path and give some output to the user
print, 'Searching for FITS files in', find_path, '...'
cd, find_path; change directories to where we're looking for FITS files
files = FILE_SEARCH('*.fits', COUNT=filecount)
print, 'Found ', filecount, ' FITS files!'

if not keyword_set(save_path) then save_path = find_path
if not keyword_set(type) then type = 'list'

; Use lists to take advantage of faster appending over arrays (we'll convert back to arrays later)
print, 'Initializing...'

; Check for the type selected and then initialize the list or array
if type eq 'list' then cube = list()
if type eq 'array' then cube = []
if type ne 'list' and type ne 'array' then begin
   print, 'Error: type needs to be "list" or "array"'
   stop
endif

; k counts what image we're on, which is used to know when to add the coadded frame to the cube
k = 1
; Start looping through each image
print, 'Initialized successfully, beginning cube creation loop...'
for ii = 0, filecount-1 do begin
   print, 'File index', ii, ' / '+string(filecount-1)
   print, 'Filename: ', files[ii], newline
   frame = readfits(files[ii])
   if type eq 'list' then cube.Add, [[frame]]
   if type eq 'array' then cube = [ [[cube]], [[frame]] ]
   
   ; print number of frames
   if type eq 'list' then print, 'Number of frames in cube: ', n_elements(cube)
   if type eq 'array' then print, 'Number of frames in cube: ', (size(cube))[3]
   k += 1
endfor; ii = start_frame for loop

; Writefits takes an array as input, so we'll need to convert our list over to an array
if type eq 'list' then begin
   print, 'Converting to array...'
   cube = cube.toArray(/TRANSPOSE, /NO_COPY)
endif

print, 'Writing cube FITS...'

;Write the cube
writefits, strcompress(save_path+'cube.fits',/rem), cube
print, 'Completed cube creation in ',(systime(/JULIAN)-starttime)*86400./60.,' minutes.'

end