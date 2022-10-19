function REDUCE_MATRIX, matrix

;Function to get rid of NaN values in 2d matrix
;Assumes that the size of finite array is rectangular i.e. that the
;NaN are the same in each image vector

n_col = (size(matrix))(1)     ; Number columns
n_row = (size(matrix))(2)     ; Number rows

index = where(finite(matrix) eq 1)

data = matrix[index]

if n_elements(data) mod n_col ne 0 then begin

   print, 'FATAL ERROR!!! IMAGE VECTORS CONTAIN ASYMMETRIC NAN CONTENT'

   return, -1

endif

index_2 = index-index[0]  ; now index starts from 0

final = fltarr(n_col,n_elements(data)/n_col)

for ii = 0, n_elements(data)-1 do begin
   
   final[index_2[ii]] = data[ii]

endfor  

return, final

end
