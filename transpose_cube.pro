pro transpose_cube, data_cube, order

compile_opt idl2

print, 'Reading in data cube...'
obj_cube = readfits(data_cube, head)
print, 'Cube read!'

print, 'Transposing data cube...'
transposed_cube = transpose(obj_cube, order)
print, 'Cube transposed!'

print, 'Writing transposed cube...'
writefits, 'transposed_cube.fits', transposed_cube
print, 'Transposed cube written!'

end
