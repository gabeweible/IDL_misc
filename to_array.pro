pro to_array
compile_opt idl2

print, 'Restoring new cube list...'
restore, '~/OneDrive/Research/HII1348/Testing/HII1348_new_cube_list.sav'
print, 'New cube list restored! Converting to array and reassigning...'
;/TRANSPOSE will give us the correct dimensions, /NO_COPY will just move over the list elements to the array, not copying them
new_cube = new_cube.toArray(/TRANSPOSE, /NO_COPY)
print, 'New cube converted to array and reassigned! Writing FITS file...'
writefits, '~/OneDrive/Research/HII1348/Testing/HII1348_cube_skysub.fits', new_cube
print, 'FITS file written!'
print, 'Finished!'

end