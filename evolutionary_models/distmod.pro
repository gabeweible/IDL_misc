function distmod,d,quiet=quiet

;takes an input (d) in parsecs and calculates magnitude modulus

md=(5. * ALOG10(d)) - 5.


return, md

end
