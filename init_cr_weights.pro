; Implementation of the Fixsen et al. (2000) CR rejection algorithm
pro init_cr_weights, bits_to_keep, gain, read_variance, sigma_cr, debug=debug
; Initialize the weight tables for CR rejection as per Fixsen et al. 2000
common cr_weights, WT, a, f, v, q
if not keyword_set(debug) then debug=0

; Constants from the paper
E = 0.0024788  ; Smallest S/N ratio bin
K = 8          ; Number of signal (S) cuts
M = 64         ; Max number of reads

; Initialize the WT array
WT = dblarr(K, M, M/2+1)

; Calculate parameters
; For LMIRCam data with 3 samples (reset + 2 reads)
N = 3  ; Number of reads

; Initialize variables used for calculations
Z = dblarr(M)
Y = dblarr(M)

; For each S/N ratio bin
for little_k = 0, K-1 do begin
    little_s = E * exp(little_k)  ; S/N ratio
    WT[little_k, 0, 0] = 0.0
    
    ; Reset variables for this iteration
    Z[*] = 0.0
    Y[*] = 0.0
    little_x = 0.0
    little_y = 0.0
    
    ; For all ramp lengths
    for little_n = 1, N-1 do begin
        W = reform(WT[little_k, little_n, *])
        
        ; New x, new last Y
        little_x = 1.0/(2.0 + little_s - little_x)
        Y[little_n] = little_x
        little_y++
        
        ; Update row and sum
        little_z = 0.0
        for i = 1, little_n-1 do begin
            Y[i] *= little_x
            Z[i] += little_y * Y[i]
            little_z += Z[i]
        endfor
        
        ; Sum weight
        little_y = little_x
        for i = 1, little_n-1 do little_y += Y[i]
        
        ; Differences - FIXED: using correct variable naming
        for i = 0, (little_n+1)/2-1 do W[i] = (Z[i] - Z[i+1]) * gain * (N-1)
        
        ; Final weight
        Z[little_n] = little_y
        W[(little_n+1)/2] = little_z + little_y
    endfor
    
    ; Renormalize last row
    little_x = 1.0/(little_z + little_y)
    for i = 0, N/2-1 do W[i] *= little_x
endfor

; Set the shared variables for use in the CR rejection routine
a = bits_to_keep * sqrt(gain)                       ; Renormalize Bits
f = 3.0 * read_variance / (N * gain * (N+1) * (N+2))  ; Output offset
q = sigma_cr / gain                                ; Renormalize cutoff
v = 2.0 * q * read_variance * (N*N + 1) / (gain * N * (N-1))  ; Renormalize variance

if debug eq 1 then begin
    print, '--- Weight Table Initialization ---'
    print, 'a value:', a
    print, 'f value:', f
    print, 'v value:', v
    print, 'q value:', q
    print, 'a * sqrt(f) =', a * sqrt(f)  ; Check if this equals 2.07954
    print, '--- End Weight Table Info ---'
endif

end