; routine to make a text file containing all the IDL-defined procedure
; and function names, for incorporation into the text module for TextWrangler
; Note: these are inbuilt routines, not those written in the IDL language

; 20060621 WJP:	initial write

PRO GET_ROUTINE_NAMES, outputFile, LCASE = lcase

IF (N_ELEMENTS(outputFile) EQ 0) THEN outputFile = DIALOG_PICKFILE(/write)

IF (N_ELEMENTS(outputFile) EQ 0) THEN BEGIN
	MESSAGE, 'No output file specified, exiting', /INFORMATIONAL
	RETURN
ENDIF

IF ~FILE_TEST(FILE_DIRNAME(outputFile), /DIRECTORY) THEN BEGIN
	MESSAGE, 'Unable to access directory: ' + file_dirname(outputFile), /INFORMATIONAL
	RETURN
ENDIF

OPENW, lun, outputFile, /GET_LUN

procedures = ROUTINE_INFO(/SYSTEM)
functions  = ROUTINE_INFO(/SYSTEM, /FUNCTIONS)

IF KEYWORD_SET(lcase) THEN procedures = STRLOWCASE(TEMPORARY(procedures))
IF KEYWORD_SET(lcase) THEN functions  = STRLOWCASE(TEMPORARY(functions))

PRINTF, lun, '<!-- System Procedures -->'
FOR ii = 0, N_ELEMENTS(procedures) - 1 do $
	PRINTF, lun, '<string>' + procedures[ii] + '</string>'

PRINTF, lun, ''
PRINTF, lun, '<!-- System Functions -->'
FOR ii = 0, N_ELEMENTS(functions) - 1 do $
	PRINTF, lun, '<string>' + functions[ii] + '</string>'
FREE_LUN, lun

;	PRINTF, lun, '<!-- Library Procedures -->'
;	PRINTF, lun, '<!-- Library Functions -->'

MESSAGE, 'System procedures and functions written to ' + outputFile, /INFORMATIONAL

RETURN

END
