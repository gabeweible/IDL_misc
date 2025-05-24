pro procyon_conversions

; Besselian dates, astrometry from:
;https://ui.adsabs.harvard.edu/abs/2015ApJ...813..106B/abstract
besselian_dates_arr = [1995.1745, 1997.1958, 1997.3747, 1997.9072, 1998.8257,$
	1999.8342, 2000.9093, 2001.8839, 2002.8537, 2003.8220, 2004.8629, 2005.9040,$
	2006.8046, 2007.8011, 2011.1040, 2012.1877, 2013.0947, 2014.7038, 2011.1040,$
	2012.1877, 2013.0947, 2014.7038, 2011.1040, 2012.1877, 2013.0947, 2014.7038]
	
besselian_ground_based = [1896.930, 1897.000, 1897.160, 1897.821, 1898.050,$
	1898.129, 1898.189, 1898.213, 1898.282, 1898.880, 1899.073, 1899.960,$
	1900.055, 1900.236, 1900.295, 1901.200, 1901.300, 1901.883, 1902.214,$
	1902.241, 1902.241, 1902.960, 1903.154, 1904.294, 1904.795, 1905.570,$
	1909.162, 1909.298, 1910.025, 1911.060, 1911.069, 1913.162, 1928.824,$
	1932.272, 1932.277, 1957.840, 1957.853, 1986.254, 1995.090]; removed those rejected by their orbital solution

; Object for Procyon B is always '1'
obj_arr = make_array((size(besselian_dates_arr))[1], /integer, value=1)
obj_arr_ground_based = make_array((size(besselian_ground_based))[1], /integer, value=1)

; nominal projected separation in arcsec
sep_arr_arcsec = [4.9389, 4.7851, 4.7651, 4.7058, 4.5973, 4.4583, 4.2809, 4.0859,$
	3.8584, 3.5988, 3.2840, 2.9293, 2.6266, 2.3452, 2.6431, 3.0130, 3.3154, 3.7986,$
	2.6381, 3.0144, 3.3166, 3.7946, 2.6406, 3.0137, 3.3160, 3.7966]
	
sep_arr_ground_based = [4.63, 4.83, 4.65, 4.66, 4.75, 4.78, 4.57, 4.83,$
	4.50, 4.97, 4.91, 4.88, 5.09, 4.83,$
	4.60, 5.13, 5.00, 5.06, 5.39, 5.34, 5.11,$
	5.33, 5.16, 4.93, 5.36,$
	5.14, 5.26, 5.04, 5.21, 4.70, 4.69, 5.09,$
	2.07,$
	3.57, 3.96, 4.554, 4.573,$
	5.10,$
	5.12]

sep_arr_mas = sep_arr_arcsec * 1000 ; convert arcseconds to mas
sep_ground_based_mas = sep_arr_ground_based * 1000 ;arcseconds to mas

; error on projected separation in arcsec
sep_err_arr_arcsec = [0.0044, 0.0047, 0.0040, 0.0030, 0.0028, 0.0027, 0.0026, 0.0032,$
	0.0029, 0.0035, 0.0032, 0.0027, 0.0027, 0.0027, 0.0047, 0.0040, 0.0040, 0.0040,$
	0.0048, 0.0040, 0.0040, 0.0040, 0.0034, 0.0028, 0.0028, 0.0028]
	
; adopted ground-based uncertainties from Bond+ 2015
sep_err_ground_based = make_array((size(sep_arr_ground_based))[1], /float, value=0.187)

sep_err_arr_mas = sep_err_arr_arcsec * 1000 ; convert arcseconds to mas
sep_err_ground_based_mas = sep_err_ground_based * 1000

; position angle in degrees
pa_arr = [42.977, 53.997, 55.022, 58.027, 63.499, 69.771, 76.977, 84.147, 91.939,$
	100.787, 112.092, 125.956, 140.997, 161.715, 240.339, 257.721, 269.417, 285.719,$
	240.205, 257.752, 269.447, 285.723, 240.272, 257.737, 269.432, 285.721]
	
pa_arr_ground_based = [320.92, 321.62, 320.32, 324.92, 325.12, 327.52, 327.12, 326.52,$
	325.52, 331.11, 331.11, 335.01, 336.54, 338.81,$
	332.91, 339.00, 338.40, 343.99, 345.40, 345.00, 347.00,$
	354.09, 351.52, 355.69, 357.87,$
	8.68, 22.97, 22.96, 26.71, 29.10, 29.05, 43.00,$
	231.06,$
	278.54, 276.04, 63.51, 64.06,$
	356.27,$
	41.00]

; error on position angle in degrees
pa_err_arr = [0.053, 0.059, 0.051, 0.039, 0.038, 0.039, 0.039, 0.049, 0.047, 0.060,$
	0.060, 0.058, 0.065, 0.074, 0.105, 0.078, 0.071, 0.062, 0.106, 0.078, 0.071,$
	0.063, 0.074, 0.055, 0.050, 0.044]
	
; Assumes equal uncertainties in delta RA and delta DEC
; This means that the fractional uncert. in SEP is the PA uncert. in rad
pa_err_ground_based_rad = sep_err_ground_based_mas / sep_ground_based_mas
pa_err_ground_based_deg = pa_err_ground_based_rad * 180/!DPI; Rad to deg 

; convert Besselian dates to MJD
mjd_arr = JBEPOCH(/B, besselian_dates_arr, /TO_DAY, /MJD)  ;; Besselian Epoch to MJD
mjd_ground_based = jbepoch(/B, besselian_ground_based, /TO_DAY, /MJD)


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
	
write_csv, '/Users/gabeweible/OneDrive/research/procyon/procyon_ground_based_mjd.csv',$
	mjd_ground_based, obj_arr_ground_based, sep_ground_based_mas,$
	sep_err_ground_based_mas, pa_arr_ground_based, pa_err_ground_based_deg,$
	header=['epoch', 'object', 'sep', 'sep_err', 'pa', 'pa_err']

end
