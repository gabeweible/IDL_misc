PRO combine_parang_files, obj, base_dir, output_dir, nod_min, nod_max
; Base_dir is where the parang.sav files are
; output_dir is the full path to output the combined results

  ; Define base directory
 ; base_dir = '/Users/gweible/OneDrive - University of Arizona/research/Alcor/macbook_5/'
 ; output_dir = '/Users/gweible/OneDrive - University of Arizona/research/Alcor/macbook_5/processed_right/'
  
  ; Initialize arrays to store combined data
  dith1_angles_combined = []
  dith1_dits_combined = []
  dith2_angles_combined = []
  dith2_dits_combined = []
  
  print, 'Starting parallactic angle combination process...'
  print, '=================================================='
  
  ; Process NOD_A files (dith1 - even nods)
  print, ''
  print, 'Processing NOD_A files (dith1 - even nods)...'
  print, '----------------------------------------------'
  
  nod_count_dith1 = 0
  for nod = nod_min, nod_max, 2 do begin  ; Even nods only
    
    ; Construct filename
    filename = string(format='("*_NOD_A_nod", I02, "_grp00_parang.sav")', nod)
    full_path = base_dir + filename
    
    ; Check if file exists
    if file_test(full_path) then begin
      print, string(format='("Reading: ", A)', filename)
      
      ; Restore the IDL save file
      restore, full_path
      
      ; Check if variables exist
      if n_elements(angles_arr) gt 0 and n_elements(dits_arr) gt 0 then begin
        print, string(format='("  - angles_arr length: ", I0)', n_elements(angles_arr))
        print, string(format='("  - dits_arr length: ", I0)', n_elements(dits_arr))
        
        ; Combine arrays
        if nod_count_dith1 eq 0 then begin
          dith1_angles_combined = angles_arr
          dith1_dits_combined = dits_arr
        endif else begin
          dith1_angles_combined = [dith1_angles_combined, angles_arr]
          dith1_dits_combined = [dith1_dits_combined, dits_arr]
        endelse
        
        nod_count_dith1++
        print, string(format='("  ✓ Successfully processed nod", I02)', nod)
        
      endif else begin
        print, string(format='("  ✗ Error: Required variables not found in nod", I02)', nod)
      endelse
      
    endif else begin
      print, string(format='("  - File not found: ", A)', filename)
    endelse
    
  endfor
  
  print, string(format='("Total NOD_A files processed: ", I0)', nod_count_dith1)
  if nod_count_dith1 gt 0 then begin
    print, string(format='("Combined dith1 angles length: ", I0)', n_elements(dith1_angles_combined))
    print, string(format='("Combined dith1 dits length: ", I0)', n_elements(dith1_dits_combined))
  endif
  
  ; Process NOD_B files (dith2 - odd nods)
  print, ''
  print, 'Processing NOD_B files (dith2 - odd nods)...'
  print, '----------------------------------------------'
  
  nod_count_dith2 = 0
  for nod = nod_min+1, nod_max, 2 do begin  ; Odd nods only
    
    ; Construct filename
   filename = string(format='("*_NOD_B_nod", I02, "_grp00_parang.sav")', nod)
    full_path = base_dir + filename
    
    ; Check if file exists
    if file_test(full_path) then begin
      print, string(format='("Reading: ", A)', filename)
      
      ; Restore the IDL save file
      restore, full_path
      
      ; Check if variables exist
      if n_elements(angles_arr) gt 0 and n_elements(dits_arr) gt 0 then begin
        print, string(format='("  - angles_arr length: ", I0)', n_elements(angles_arr))
        print, string(format='("  - dits_arr length: ", I0)', n_elements(dits_arr))
        
        ; Combine arrays
        if nod_count_dith2 eq 0 then begin
          dith2_angles_combined = angles_arr
          dith2_dits_combined = dits_arr
        endif else begin
          dith2_angles_combined = [dith2_angles_combined, angles_arr]
          dith2_dits_combined = [dith2_dits_combined, dits_arr]
        endelse
        
        nod_count_dith2++
        print, string(format='("  ✓ Successfully processed nod", I02)', nod)
        
      endif else begin
        print, string(format='("  ✗ Error: Required variables not found in nod", I02)', nod)
      endelse
      
    endif else begin
      print, string(format='("  - File not found: ", A)', filename)
    endelse
    
  endfor
  
  print, string(format='("Total NOD_B files processed: ", I0)', nod_count_dith2)
  if nod_count_dith2 gt 0 then begin
    print, string(format='("Combined dith2 angles length: ", I0)', n_elements(dith2_angles_combined))
    print, string(format='("Combined dith2 dits length: ", I0)', n_elements(dith2_dits_combined))
  endif
  
  ; Save combined data for dith1
  if nod_count_dith1 gt 0 then begin
    print, ''
    print, 'Saving dith1 combined data...'
    
    ; Assign to output variable names
    dx_dith1_angles_array = dith1_angles_combined
    dx_dith1_dits_array = dith1_dits_combined
    
    ; Create output directory if it doesn't exist
    dith1_dir = output_dir + 'dith1/'
    if ~file_test(dith1_dir, /directory) then file_mkdir, dith1_dir
    
    ; Save to file
    dith1_output_file = dith1_dir + obj + '_parang.sav'
    save, dx_dith1_angles_array, dx_dith1_dits_array, filename=dith1_output_file
    writefits, dith1_dir + obj+'_angles.fits', dx_dith1_angles_array
    writefits, dith1_dir + obj+'_dits.fits', dx_dith1_dits_array
    print, string(format='("✓ Saved dith1 data to: ", A)', dith1_output_file)
    print, string(format='("  - dx_dith1_angles_array length: ", I0)', n_elements(dx_dith1_angles_array))
    print, string(format='("  - dx_dith1_dits_array length: ", I0)', n_elements(dx_dith1_dits_array))
  endif
  
  ; Save combined data for dith2
  if nod_count_dith2 gt 0 then begin
    print, 'Saving dith2 combined data...'
    
    ; Assign to output variable names
    dx_dith2_angles_array = dith2_angles_combined
    dx_dith2_dits_array = dith2_dits_combined
    
    ; Create output directory if it doesn't exist
    dith2_dir = output_dir + 'dith2/'
    if ~file_test(dith2_dir, /directory) then file_mkdir, dith2_dir
    
    ; Save to file
    dith2_output_file = dith2_dir + obj+'_parang.sav'
    save, dx_dith2_angles_array, dx_dith2_dits_array, filename=dith2_output_file
    writefits, dith2_dir + obj+'_angles.fits', dx_dith2_angles_array
    writefits, dith2_dir + obj+'_dits.fits', dx_dith2_dits_array
    print, string(format='("✓ Saved dith2 data to: ", A)', dith2_output_file)
    print, string(format='("  - dx_dith2_angles_array length: ", I0)', n_elements(dx_dith2_angles_array))
    print, string(format='("  - dx_dith2_dits_array length: ", I0)', n_elements(dx_dith2_dits_array))
  endif
  
  print, ''
  print, '=================================================='
  print, 'Parallactic angle combination process completed!'
  
  ; Print summary
  print, ''
  print, 'SUMMARY:'
  print, string(format='("- NOD_A files processed: ", I0)', nod_count_dith1)
  print, string(format='("- NOD_B files processed: ", I0)', nod_count_dith2)
  if nod_count_dith1 gt 0 then begin
    print, string(format='("- dith1 total angles: ", I0)', n_elements(dx_dith1_angles_array))
    print, string(format='("- dith1 total dits: ", I0)', n_elements(dx_dith1_dits_array))
  endif
  if nod_count_dith2 gt 0 then begin
    print, string(format='("- dith2 total angles: ", I0)', n_elements(dx_dith2_angles_array))
    print, string(format='("- dith2 total dits: ", I0)', n_elements(dx_dith2_dits_array))
  endif

END