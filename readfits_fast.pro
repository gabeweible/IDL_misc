function READFITS_fast, filename, header, heap, CHECKSUM=checksum, $ 
                   COMPRESS = compress, HBUFFER=hbuf, EXTEN_NO = exten_no, $
                   NOSCALE = noscale, NSLICE = nslice, $
                   NO_UNSIGNED = no_unsigned,  NUMROW = numrow, $
                   POINTLUN = pointlun, SILENT = silent, STARTROW = startrow, $
                   NaNvalue = NaNvalue, FPACK = fpack, UNIXpipe=unixpipe

  compile_opt idl2
  On_IOerror, BAD

; Check for filename input
   if N_params() LT 1 then begin                
      print,'Syntax - im = READFITS( filename, [ h, heap, /NOSCALE, /SILENT,'
      print,'                 EXTEN_NO =, STARTROW = , NUMROW=, NSLICE = ,'
      print,'                 HBUFFER = ,/NO_UNSIGNED, /CHECKSUM, /COMPRESS]'
      return, -1
   endif
   
Catch, theError
if theError NE 0 then begin
    Catch,/Cancel
    void = cgErrorMsg(/quiet)
    return, -1
endif

; Set default keyword values upfront to avoid repeated checks
unitsupplied = size(filename,/TNAME) EQ 'STRING' ? 0 : 1
silent = keyword_set(SILENT)
do_checksum = keyword_set(CHECKSUM)
exten_no = N_elements(exten_no) EQ 0 ? 0 : exten_no
unixpipe = N_elements(unixpipe) EQ 0 ? 0 : unixpipe

; Pre-define variables that will be used in loops
header_block = 0L
i = 0L

; Set up file handling - optimized by reducing string operations and combining conditions
if unitsupplied then begin
    unit = filename
endif else begin
    len = strlen(filename)
    if len EQ 0 then begin
        filename = dialog_pickfile(filter=['*.fit*;*.fts*;*.img*'], $
                   title='Please select a FITS file',/must_exist)
        len = strlen(filename)
    endif
    
    ; Determine file type with fewer string operations
    ext = strlowcase(strmid(filename, len-3, 3))
    gzip = (ext EQ '.gz') || (ext EQ 'ftz')
    compress = keyword_set(compress) || gzip
    unixZ = (strmid(filename, len-2, 2) EQ '.Z')
    bzip = (ext EQ 'bz2')
    fcompress = keyword_set(fpack) || (ext EQ '.fz')
    unixpipe = unixZ || fcompress || bzip
    
    ; Open file - combined error handling
    openr, unit, filename, ERROR=error, /get_lun, $
           COMPRESS=compress, /swap_if_little_endian
    
    if error NE 0 then begin
        message, /con, ' ERROR - Unable to locate file ' + filename
        return, -1
    endif
    
    ; Handle compressed files more efficiently
    if unixZ || bzip then begin
        free_lun, unit
        cmd = bzip ? 'bunzip2' : 'gzip'
        spawn, cmd + ' -cd ' + filename, unit=unit
    endif else if fcompress then begin
        free_lun, unit
        spawn, 'funpack -S ' + filename, unit=unit, /sh
        if eof(unit) then begin
            message, 'Error spawning FPACK decompression', /CON
            free_lun, unit
            return, -1
        endif
    endif
endelse

; Handle POINTLUN once
if N_elements(POINTLUN) GT 0 then mrd_skip, unit, pointlun

; Optimize header buffer size calculation
if N_elements(hbuf) EQ 0 then begin
    hbuf = 180
endif else begin
    remain = hbuf mod 36
    if remain GT 0 then hbuf += (36 - remain)
endelse

; Improve extension loop by preallocating memory and reducing operations
for ext = 0L, exten_no do begin
    ; Pre-allocate header array to avoid resizing
    if (ext EQ exten_no) then begin
        header = strarr(hbuf)
    endif else begin
        header = strarr(36)
    endelse
    
    header_block = 0L
    i = 0L
    found_end = 0
    
    ; Read header more efficiently
    while ~found_end do begin
        if EOF(unit) then begin
            message, /CON, 'EOF encountered attempting to read extension ' + strtrim(ext, 2)
            if ~unitsupplied then free_lun, unit
            return, -1
        endif
        
        ; Read a block of header data (80 chars × 36 lines)
        block = string(replicate(32b, 80, 36))
        readu, unit, block
        header_block++
        
        ; Check for invalid characters - vectorized operation
        w = where(strlen(block) NE 80, Nbad)
        if (Nbad GT 0) then begin
            if ~silent then message, 'Warning-Invalid characters in header', /INF
            block[w] = string(replicate(32b, 80))
        endif
        
        ; Find END marker - using single WHERE operation
        w = where(strcmp(block, 'END     ', 8), Nend)
        found_end = (Nend GT 0)
        
        ; Store header efficiently
        if (header_block EQ 1) || (ext EQ exten_no) then begin
            if found_end then begin
                if header_block EQ 1 then begin
                    header = block[0:w[0]]
                endif else begin
                    if i+w[0]+1 LE n_elements(header) then begin
                        header[i:i+w[0]] = block[0:w[0]]
                    endif else begin
                        temp = header
                        header = [temporary(temp), block[0:w[0]]]
                    endelse
                endelse
            endif else begin
                if i+36 LE n_elements(header) then begin
                    header[i:i+35] = block
                endif else begin
                    temp = header
                    header = [temporary(temp), block]
                endelse
                i += 36
            endelse
        endif
        
        ; Check SIMPLE keyword only when needed (first extension, no pointlun)
        if (ext EQ 0) && ~((N_elements(pointlun) GT 0) || unitsupplied) then begin
            if ~found_end && (header_block EQ 1) && (strmid(header[0], 0, 8) NE 'SIMPLE  ') then begin
                message, /CON, 'ERROR - Header does not contain required SIMPLE keyword'
                if ~unitsupplied then free_lun, unit
                return, -1
            endif
        endif
    endwhile
    
    ; Get parameters more efficiently
    bitpix = sxpar(header, 'BITPIX')
    byte_elem = abs(bitpix) / 8
    naxis = sxpar(header, 'NAXIS')
    gcount = sxpar(header, 'GCOUNT') > 1
    pcount = sxpar(header, 'PCOUNT')
    
    ; Calculate array dimensions more efficiently
    if naxis GT 0 then begin
        dims = sxpar(header, 'NAXIS*')
        ndata = product(dims[0:naxis-1], /integer)
    endif else begin
        ndata = 0L
    endelse
    
    nbytes = byte_elem * gcount * (pcount + ndata)
    
    ; Skip to next extension efficiently
    if ext LT exten_no then begin
        nrec = long64((nbytes + 2879) / 2880)
        if nrec GT 0 then mrd_skip, unit, nrec * 2880L
    endif
endfor

; Set up data type more efficiently with direct mapping
case BITPIX of
    8:   IDL_type = 1          ; Byte
    16:  IDL_type = 2          ; 16 bit integer
    32:  IDL_type = 3          ; 32 bit integer
    64:  IDL_type = 14         ; 64 bit integer
    -32: IDL_type = 4          ; Float
    -64: IDL_type = 5          ; Double
    else: begin
        message, /CON, 'ERROR - Illegal value of BITPIX (= ' + strtrim(bitpix, 2) + ') in FITS header'
        if ~unitsupplied then free_lun, unit
        return, -1
    end
endcase

; Handle zero data case more efficiently
if nbytes EQ 0 then begin
    if ~silent then message, "FITS header has NAXIS or NAXISi = 0, no data array read", /CON
    if do_checksum then begin
        result = FITS_TEST_CHECKSUM(header, data, ERRMSG=errmsg)
        if ~silent then begin
            case result of
                1: message, /INF, 'CHECKSUM keyword in header is verified'
                -1: message, /CON, errmsg
                else:
            endcase
        endif
    endif
    if ~unitsupplied then free_lun, unit
    return, -1
endif

; Check for FITS extensions efficiently
groups = sxpar(header, 'GROUPS')
if groups && ~silent then message, 'WARNING - FITS file contains random GROUPS', /INF

; Handle row selection efficiently
startrow = N_elements(STARTROW) EQ 0 ? 0 : startrow
if naxis GE 2 then nrow = dims[1] else nrow = ndata
numrow = N_elements(NUMROW) EQ 0 ? nrow : numrow

; Optimize checksum handling
if do_checksum && ((startrow GT 0) || (numrow LT nrow) || (N_elements(nslice) GT 0)) then begin
    if ~silent then message, /CON, 'Warning - CHECKSUM not applied when STARTROW, NUMROW or NSLICE is set'
    do_checksum = 0
endif

; Extension handling
if exten_no GT 0 then begin
    xtension = strtrim(sxpar(header, 'XTENSION', Count=N_ext), 2)
    if (N_ext EQ 0) && ~silent then message, /INF, 'WARNING - Header missing XTENSION keyword'
endif

; Optimize row and slice selection
if (startrow NE 0) || (numrow NE nrow) then begin
    if startrow GE dims[1] then begin
        message, 'ERROR - Specified starting row ' + strtrim(startrow, 2) + $
                ' but only ' + strtrim(dims[1], 2) + ' rows in extension', /CON
        if ~unitsupplied then free_lun, unit
        return, -1
    endif
    
    dims[1] = dims[1] - startrow
    dims[1] = dims[1] < numrow
    sxaddpar, header, 'NAXIS2', dims[1]
    
    if startrow GT 0 then mrd_skip, unit, byte_elem * startrow * dims[0]
endif else if N_elements(NSLICE) EQ 1 then begin
    ; Optimize slice handling with fewer operations
    ldim = naxis - 1
    lastdim = dims[ldim]
    
    while lastdim EQ 1 && ldim GT 0 do begin
        ldim--
        lastdim = dims[ldim]
    endwhile
    
    if nslice GE lastdim then begin
        message, /CON, 'ERROR - Value of NSLICE must be less than ' + strtrim(lastdim, 2)
        if ~unitsupplied then free_lun, unit
        return, -1
    endif
    
    ; Update dimensions more efficiently
    dims = dims[0:ldim-1]
    for i = ldim, naxis-1 do sxdelpar, header, 'NAXIS' + strtrim(i+1, 2)
    naxis = ldim
    sxaddpar, header, 'NAXIS' + strtrim(ldim, 2), 1
    
    ndata = ndata / lastdim
    nskip = long64(nslice) * ndata * byte_elem
    if ndata GT 0 then mrd_skip, unit, nskip
endif

; Print array info efficiently
if ~silent then begin
    if exten_no GT 0 then message, 'Reading FITS extension of type ' + xtension, /INF
    
    if N_elements(dims) EQ 1 then begin
        st = 'Now reading ' + strtrim(dims, 2) + ' element vector'
    endif else begin
        st = 'Now reading ' + strjoin(strtrim(dims, 2), ' by ') + ' array'
    endelse
    
    if (exten_no GT 0) && (pcount GT 0) then st += ' + heap area'
    message, /INF, st
endif

; Read data efficiently - using MAKE_ARRAY with /NOZERO for speed
data = make_array(DIM=dims, TYPE=IDL_type, /NOZERO)
readu, unit, data

; Handle byte swapping efficiently
if unixpipe then swap_endian_inplace, data, /swap_if_little

; Heap handling optimized
if (exten_no GT 0) && (pcount GT 0) then begin
    theap = sxpar(header, 'THEAP')
    skip = theap - N_elements(data)
    
    if skip GT 0 then begin
        temp = bytarr(skip, /nozero)
        readu, unit, skip
    endif
    
    heap = bytarr(pcount * gcount * byte_elem)
    readu, unit, heap
    
    if do_checksum then begin
        result = fits_test_checksum(header, [temporary(data), heap], ERRMSG=errmsg)
    endif
endif else if do_checksum then begin
    result = fits_test_checksum(header, data, ERRMSG=errmsg)
endif

; Close file if we opened it
if ~unitsupplied then free_lun, unit

; Checksum reporting
if do_checksum && ~silent then begin
    case result of
        1: message, /INF, 'CHECKSUM keyword in header is verified'
        -1: message, /CON, 'CHECKSUM ERROR! ' + errmsg
        else:
    endcase
endif

; Optimize scaling operations
do_scale = ~keyword_set(NOSCALE)
if (do_scale && (exten_no GT 0)) then do_scale = xtension EQ 'IMAGE'

if do_scale then begin
    ; Get scaling parameters once
    if bitpix GT 0 then begin
        blank = sxpar(header, 'BLANK', Count=N_blank)
    endif else begin
        N_blank = 0
    endelse
    
    Bscale = sxpar(header, 'BSCALE', Count=N_bscale)
    Bzero = sxpar(header, 'BZERO', Count=N_Bzero)
    
    if (N_blank GT 0) && ((N_bscale GT 0) || (N_Bzero GT 0)) then begin
        sxaddpar, header, 'O_BLANK', blank, ' Original BLANK value'
    endif
    
    ; Handle unsigned integer efficiently with combined conditions
    unsgn = 0
    if ~keyword_set(No_Unsigned) then begin
        no_bscale = (Bscale EQ 1) || (N_bscale EQ 0)
        unsgn_int = (bitpix EQ 16) && (Bzero EQ 32768) && no_bscale
        unsgn_lng = (bitpix EQ 32) && (Bzero EQ 2147483648) && no_bscale
        unsgn_lng64 = (bitpix EQ 64) && (Bzero EQ 9223372036854775808) && no_bscale
        
        unsgn = unsgn_int || unsgn_lng || unsgn_lng64
    endif
    
    ; Perform scaling operations more efficiently
    if unsgn then begin
        if unsgn_int then begin
            data = uint(temporary(data)) - 32768US
            if N_blank then blank = uint(blank) - 32768US
        endif else if unsgn_lng then begin
            data = ulong(temporary(data)) - 2147483648UL
            if N_blank then blank = ulong(blank) - 2147483648UL
        endif else begin
            offset = ulong64(9223372036854775808)
            data = ulong64(temporary(data)) - offset
            if N_blank then blank = ulong64(blank) - offset
        endelse
        
        if N_blank then sxaddpar, header, 'BLANK', blank
        sxaddpar, header, 'BZERO', 0
        sxaddpar, header, 'O_BZERO', Bzero, ' Original BZERO Value'
    endif else begin
        ; Apply BSCALE and BZERO more efficiently
        scale_applied = 0
        
        if (N_Bscale GT 0) && (Bscale NE 1.0) then begin
            ; Apply scaling with minimal type conversions
            if size(Bscale, /TNAME) NE 'DOUBLE' then begin
                data = temporary(data) * float(Bscale)
            endif else begin
                data = temporary(data) * Bscale
            endelse
            
            if N_blank then blank *= Bscale
            sxaddpar, header, 'BSCALE', 1.0
            sxaddpar, header, 'O_BSCALE', Bscale, ' Original BSCALE Value'
            scale_applied = 1
        endif
        
        if (N_Bzero GT 0) && (Bzero NE 0) then begin
            ; Apply zero point with minimal type conversions
            if size(Bzero, /TNAME) NE 'DOUBLE' then begin
                data = temporary(data) + float(Bzero)
            endif else begin
                data = temporary(data) + Bzero
            endelse
            
            if N_blank then blank += Bzero
            sxaddpar, header, 'BZERO', 0.0
            sxaddpar, header, 'O_BZERO', Bzero, ' Original BZERO Value'
            scale_applied = 1
        endif
    endelse
    
    ; Update BLANK value if needed
    if N_blank then sxaddpar, header, 'BLANK', blank
endif

; Handle NaN values efficiently
if n_elements(nanvalue) EQ 1 then begin
    w = where(finite(data, /nan), count)
    if count GT 0 then data[w] = nanvalue
endif

; Return data efficiently
return, data

; Error handler
BAD:
    print, !ERROR_STATE.MSG
    if (~unitsupplied) && (N_elements(unit) GT 0) then free_lun, unit
    if N_elements(data) GT 0 then return, data else return, -1
end