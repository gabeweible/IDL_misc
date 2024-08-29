pro procyon_conversions

; Besselian dates, astrometry from:
;https://ui.adsabs.harvard.edu/abs/2015ApJ...813..106B/abstract
besselian_dates_arr = [1995.1745, 1997.1958, 1997.3747, 1997.9072, 1998.8257,$
	1999.8342, 2000.9093, 2001.8839, 2002.8537, 2003.8220, 2004.8629, 2005.9040,$
	2006.8046, 2007.8011, 2011.1040, 2012.1877, 2013.0947, 2014.7038, 2011.1040,$
	2012.1877, 2013.0947, 2014.7038, 2011.1040, 2012.1877, 2013.0947, 2014.7038]

; Object for Procyon B is always '1'
obj_arr = make_array((size(besselian_dates_arr))[1], /integer, value=1)

; nominal projected separation in arcsec
sep_arr_arcsec = [4.9389, 4.7851, 4.7651, 4.7058, 4.5973, 4.4583, 4.2809, 4.0859,$
	3.8584, 3.5988, 3.2840, 2.9293, 2.6266, 2.3452, 2.6431, 3.0130, 3.3154, 3.7986,$
	2.6381, 3.0144, 3.3166, 3.7946, 2.6406, 3.0137, 3.3160, 3.7966]

sep_arr_mas = sep_arr_arcsec * 1000 ; convert arcseconds to mas

; error on projected separation in arcsec
sep_err_arr_arcsec = [0.0044, 0.0047, 0.0040, 0.0030, 0.0028, 0.0027, 0.0026, 0.0032,$
	0.0029, 0.0035, 0.0032, 0.0027, 0.0027, 0.0027, 0.0047, 0.0040, 0.0040, 0.0040,$
	0.0048, 0.0040, 0.0040, 0.0040, 0.0034, 0.0028, 0.0028, 0.0028]
	
sep_err_arr_mas = sep_err_arr_arcsec * 1000 ; convert arcseconds to mas

; position angle in degrees
pa_arr = [42.977, 53.997, 55.022, 58.027, 63.499, 69.771, 76.977, 84.147, 91.939,$
	100.787, 112.092, 125.956, 140.997, 161.715, 240.339, 257.721, 269.417, 285.719,$
	240.205, 257.752, 269.447, 285.723, 240.272, 257.737, 269.432, 285.721]

; error on position angle in degrees
pa_err_arr = [0.053, 0.059, 0.051, 0.039, 0.038, 0.039, 0.039, 0.049, 0.047, 0.060,$
	0.060, 0.058, 0.065, 0.074, 0.105, 0.078, 0.071, 0.062, 0.106, 0.078, 0.071,$
	0.063, 0.074, 0.055, 0.050, 0.044]

; convert Besselian dates to MJD
mjd_arr = JBEPOCH(/B, besselian_dates_arr, /TO_DAY, /MJD)  ;; Besselian Epoch to MJD


;print, mjd_arr
;print, obj_arr
;print, sep_arr_mas
;print, sep_err_arr_mas
;print, pa_arr
;print, pa_err_arr
; output to a CSV that we can add some empty columns to for reading into orbitize!
write_csv, '/Users/gabeweible/OneDrive/research/procyon/procyon_obs_mjd.csv',$
	mjd_arr, obj_arr, sep_arr_mas, sep_err_arr_mas, pa_arr, pa_err_arr,$
	header=['epoch', 'object', 'sep', 'sep_err', 'pa', 'pa_err']

end
