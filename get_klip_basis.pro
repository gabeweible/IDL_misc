FUNCTION GET_KLIP_BASIS, image_arr, k_klip

k_klipi=k_klip ;save initial k_klip state

;Function that returns the KL basis of input matrix
;IMAGE_ARR - Data matrix : N_obj (columns) X M_attrib (rows)
;K_KLIP - Integer number of KL Basis Vectors to keep

k_klip = k_klip < (size(image_arr))[1]  ; AD: k_klip cannot be larger the number of frames available

;Data Matix "R"
R = image_arr

;Get Dimensions
R_Dimen = (size(image_arr))(0)

n_col = (size(R))(1)
n_row = (size(R))(2)

;Make sure we read the correct dimensions!
if R_Dimen eq 1 then begin

n_col = (size(R))(1)  ; size now returns col as arg. (1)
n_row = 1.0

endif

; Constants
TOLERANCE = 5.0E-3              ; Close enough to zero?

;Subtract the mean of each object from itself
for ii = 0,n_col-1 do begin 

   R[ii,*] = R[ii,*] - total(R[ii,*])/n_row

endfor

;We want the eigendecomposition of RR'
A = R # transpose(R)

eigenval = eigenql(A, EIGENVECTORS = eigenvect, /double)

n_eig = n_elements(eigenval)

eig_val_arr = fltarr(n_eig,n_eig)

for jj = 0, n_elements(eigenval)-1 do begin

   if abs(eigenval[jj]) lt TOLERANCE then eigenval[jj] = 0.0

   eig_val_arr[jj,jj] = eigenval[jj]

   if sqrt(eigenval[jj]) ne 0 then begin

      eigenvect[*,jj] = eigenvect[*,jj]/sqrt(abs(eigenval[jj])) ; normalize

   endif

endfor

new_basis = matrix_multiply(eigenvect,R,/atranspose) ; kth column is kth principle component of R

;counter = 0
;percent = 0.0

;while percent lt max_percent do begin

;   percent = percent+eigenval[counter]/total(eigenval)*100

;   counter = counter+1
   
;endwhile

;k_klip = counter

n_rows = (size(new_basis))(2)

;check if 1d representation
new_basis_Dimen = (size(new_basis))(0)
if new_basis_Dimen eq 1 then n_rows = 1


kl_basis = fltarr(k_klip,n_rows)
kl_basis[*] = !values.f_nan

for jj = 0,k_klip-1 do begin 
   
   kl_basis[jj,*] = new_basis[jj,*]
   
endfor

k_klip=k_klipi	;restore initial k_klip state

return, kl_basis

END
