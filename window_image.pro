FUNCTION WINDOW_IMAGE, image, radius, buffer_width, center_x=center_x, center_y=center_y, new_size=new_size
  ;; Purpose: Windows an image to smoothly fade to mean value beyond specified radius,
  ;;          with options to specify center and resize/crop
  ;; Inputs:
  ;;   image - 2D array containing the image
  ;;   radius - radius beyond which to begin windowing
  ;;   buffer_width - width of the transition buffer zone
  ;; Optional inputs:
  ;;   center_x - x-coordinate of center (default: image center)
  ;;   center_y - y-coordinate of center (default: image center)
  ;;   new_size - 2-element array [new_nx, new_ny] for output size (default: original size)
  ;; Returns:
  ;;   Windowed image with smooth transition to mean value at edges, centered and sized as specified
  
  ;; Get image dimensions
  sz = SIZE(image)
  nx = sz[1]
  ny = sz[2]
  
  ;; Set defaults for optional parameters
  IF N_ELEMENTS(center_x) EQ 0 THEN center_x = (nx-1)/2.0
  IF N_ELEMENTS(center_y) EQ 0 THEN center_y = (ny-1)/2.0
  IF N_ELEMENTS(new_size) EQ 0 THEN new_size = [nx, ny]
  
  new_nx = new_size[0]
  new_ny = new_size[1]
  
  ;; Create a resized/cropped version of the image centered at the specified point
  resized = FLTARR(new_nx, new_ny)
  
  ;; Calculate corners of the new image in the original image coordinates
  x_min = center_x - (new_nx-1)/2.0
  x_max = center_x + (new_nx-1)/2.0
  y_min = center_y - (new_ny-1)/2.0
  y_max = center_y + (new_ny-1)/2.0
  
  ;; Create coordinate arrays for interpolation
  x_orig = FINDGEN(nx)
  y_orig = FINDGEN(ny)
  
  ;; Create grid of coordinates for the new image
  x_new = FINDGEN(new_nx) + x_min
  y_new = FINDGEN(new_ny) + y_min
  
  ;; Use bilinear interpolation to resample the image
  resized = INTERPOLATE(image, REBIN(x_new, new_nx, new_ny) - x_min + (new_nx-1)/2.0, $
                         REBIN(REFORM(y_new, 1, new_ny), new_nx, new_ny) - y_min + (new_ny-1)/2.0, $
                         /GRID, MISSING=0.0)
  
  ;; Create coordinate grid centered on new image
  x_center = (new_nx-1)/2.0
  y_center = (new_ny-1)/2.0
  x = FINDGEN(new_nx) - x_center
  y = FINDGEN(new_ny) - y_center
  
  xx = REBIN(x, new_nx, new_ny)
  yy = REBIN(REFORM(y, 1, new_ny), new_nx, new_ny)
  r = SQRT(xx^2 + yy^2)
  
  ;; Calculate mean within the specified radius
  mask = r LE radius
  IF TOTAL(mask) GT 0 THEN BEGIN
    img_mean = TOTAL(resized * mask) / TOTAL(mask)
  ENDIF ELSE BEGIN
    img_mean = MEAN(resized)
  ENDELSE
  
  ;; Create weighting function for smooth transition
  weight = FLTARR(new_nx, new_ny) + 1.0  ; start with all 1's
  
  ;; Apply cosine taper in buffer zone
  buffer_idx = WHERE((r GT radius) AND (r LE (radius + buffer_width)), count)
  IF count GT 0 THEN BEGIN
    buffer_r = r[buffer_idx] - radius
    weight[buffer_idx] = 0.5 * (1 + COS(!PI * buffer_r / buffer_width))
  ENDIF
  
  ;; Set weight to zero beyond radius + buffer
  outside_idx = WHERE(r GT (radius + buffer_width), count)
  IF count GT 0 THEN weight[outside_idx] = 0.0
  
  ;; Create windowed image: original * weight + mean * (1-weight)
  windowed = resized * weight + img_mean * (1.0 - weight)
  
  RETURN, windowed
END