pro crop_image, image, half_width_x, half_width_y, center_px, output_file
;+
; NAME:
;       crop_iamge

; PURPOSE:
;       crop an image to width_x by width_y, centered on centerpx

; CALLING SEQUENCE:
;       crop_image, image, width_x, width_y, center_px output_file

; INPUTS:
;       image = a single 2D FITS image
; 
;       width_x = length of pixels to crop image to in the x-direction

;		  width_y = length of pixels to crop image to in the y-direction

;		  center_px = IDL array coordinate of the pixel to center the cropping on

;       output_file = string output FITS to write the masked image to
;       
; RESTRICTIONS:
;       (1) image must be a string path to the FITS image, width_x must be an
;           integer, width_y must be an integer, center_px must be an integer,
;				and output_file must be a string .fits file.
;
; EXAMPLE:
;       Crop an image to 10 pixels by 12 pixels, centered on the image index [102, 80]:

;       CROP_IMAGE, '/path/to/image.fits', 10, 12, [102, 80], '/path/to/masked/image.fits'
;
; EXTERNAL PROCEDURES USED:
;       READFITS, WRITEFITS
;
; MODIFICATION HISTORY:
;       WRITTEN, 2024 Gabriel Weible
;-

; Strict arrays (square brackets only) and 32-bit integers as default, instead of 16
compile_opt IDL2
newline = string(10B)
; Get the current time in a Julian date format
; (number of days since January 1, 4713 BC, with decimals)
start_time = systime(/JULIAN)

; read-in the unmasked image
image = readfits(image)

; crop the image to width_x pixels by width_y pixels, centered at image[center_px]
cropped_image = image[center_px[0] - half_width_x : $
	center_px[0] + half_width_x,center_px[1]-half_width_y : center_px[1]+half_width_y]
	
; write the masked image
writefits, output_file, cropped_image

end