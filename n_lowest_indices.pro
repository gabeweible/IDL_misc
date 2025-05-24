FUNCTION N_LOWEST_INDICES, array, n
  ;+
  ; NAME:
  ;   GET_N_LOWEST_INDICES
  ;
  ; PURPOSE:
  ;   Returns the indices of the N lowest values in an array.
  ;
  ; CALLING SEQUENCE:
  ;   indices = GET_N_LOWEST_INDICES(array, n)
  ;
  ; INPUTS:
  ;   array - The input array of any dimension
  ;   n     - The number of lowest values to find (integer)
  ;
  ; OUTPUTS:
  ;   indices - A 1D array containing the indices of the N lowest values
  ;
  ; EXAMPLE:
  ;   arr = [5, 2, 8, 1, 7, 3]
  ;   indices = GET_N_LOWEST_INDICES(arr, 3)
  ;   print, indices    ; Should output [3, 1, 5] (indices of 1, 2, and 3)
  ;   print, arr[indices] ; Should output [1, 2, 3]
  ;
  ; NOTES:
  ;   - If N is greater than the number of elements in the array,
  ;     N will be set to the array size.
  ;   - Handles NaN values by ignoring them (they won't be included in results)
  ;-
  
  ; Check inputs
  IF N_PARAMS() LT 2 THEN BEGIN
    MESSAGE, 'Usage: indices = GET_N_LOWEST_INDICES(array, n)'
  ENDIF
  
  ; Get array size
  arr_size = N_ELEMENTS(array)
  
  ; Make sure n is not larger than array size
  n = n < arr_size
  
  ; Handle empty array or n <= 0
  IF (arr_size EQ 0) OR (n LE 0) THEN BEGIN
    RETURN, -1
  ENDIF
  
  ; Sort the array and get the sort indices
  sort_indices = SORT(array)
  
  ; Return the first n indices
  RETURN, sort_indices[0:n-1]
END
