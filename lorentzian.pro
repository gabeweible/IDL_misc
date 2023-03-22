; NAME:
;       Lorentz_deriv
; PURPOSE:
;       Evaluate partial derivatives of 1-D generalized Lorentzian function.
;       (also known as Cauchy function).  Called by function Lorentian in loop.
; CALLING:
;       pd = Lorentz_deriv( xm, radius,psc,power, alpha,beta,gama, a_b,Logalf )
; INPUTS:
;       xi = array, indepented variable
;       parms = array, generalized Lorentzian (Cauchy) function parameters (5):
; OUTPUT:
;       Function returns matrix of partial derivatives respect to 3 parameters:
;               centroid, radius, pscale.
; HISTORY:
;       Written: Frank Varosi NASA/GSFC 1994.

function Lorentz_deriv, xm, radius, psc, power, alpha, beta, gama, a_b, Logalf

        pderiv = fltarr( N_elements( alpha ), 3 )

        w = where( alpha GT 0, np )
        a_bm2 = a_b
        if (np GT 0) then a_bm2(w) = a_b(w) / alpha(w)^2

        pderiv(0,0) = xm * a_bm2 * beta / radius^2
        w = where( gama GT 0, np )

        if (np GT 0) then pderiv(w,0) = pderiv(w,0) + $
                (power * psc^2) * ( xm(w) * Logalf(w) * a_b(w) / gama(w) )

        xm2 = xm^2
        pderiv(0,1) = beta * a_bm2 * xm2 / radius^3

        if (np GT 0) then pderiv(w,2) = -(power * psc) * $
                                ( Logalf(w) * a_b(w) * xm2(w) / gama(w) )
return, pderiv
end

;+
; NAME:
;       Lorentzian
; PURPOSE:
;       Evaluate the generalized Lorentzian function.
;                       (also known as Cauchy function).
; CALLING EXAMPLES:
;       y = Lorentzian( xi, parms )
;       y = Lorentzian( xi, parms, pderiv )
; INPUTS:
;       xi = array, indepented variable
;       parms = array, generalized Lorentzian (Cauchy) function parameters (5):
; OUTPUT (optional):
;       pderiv = matrix of partial derivatives respect to parameters.
; PROCEDURE:
;       Case of dimension breakdown.
; HISTORY:
;       Written: Frank Varosi NASA/GSFC 1992.
;       F.V.1994, added optional sixth parameter = constant offset.
;       F.V.1994, generalized to 2-D and 3-D.
;-

function Lorentzian, xi, parameters, pderiv

        szx = size( xi )
        szp = size( parameters )

        if (szp(0) EQ 2) AND (szp(2) GE 5) then begin
                parms = float( transpose( parameters ) )
                szp = size( parms )
          endif else parms = float( parameters )

        Nparm = szp(1)

        if (Nparm LT 5) then begin
                message,"need at least 5 parameters",/INFO
                return,(-1)
           endif

        if (szp(0) GT 1) OR (szx(0) GT 1) then begin
           if (szp(0) NE 2) OR (szx(0) NE 2) OR (szp(2) NE szx(2)) then begin
                message,"dimensionality of independent variable and " + $
                        "parameters do not agree",/INFO
                print,szx,szp
                return,(-1)
              endif
            ndim = szp(2)
          endif else ndim = 1

CASE ndim OF

1: BEGIN
        factor = parms(0)
        cntrd = parms(1)
        radius = parms(2)
        pscale = parms(3)
        power = parms(4)

        xm = xi - cntrd
        xa = abs( xm )
        alpha = xa/radius
        gama = xa * abs( pscale )
        beta = power * ( 1 + gama )
        a_b = ( alpha^beta ) < 1.e30
        Lcf = 1/(1 + a_b)

        if N_params() GE 3 then begin

                pderiv = fltarr( N_elements( alpha ), Nparm )
                pderiv(0,0) = Lcf
                Logalf = aLog( alpha > 1.e-33 )

                pderiv(0,1) = Lorentz_deriv( xm, radius, pscale, power, $
                                             alpha, beta, gama, a_b, Logalf )

                pderiv(0,4) = -a_b * Logalf * (1 + gama)

                df = factor * Lcf^2
                for i=1,4 do pderiv(0,i) = df * pderiv(*,i)
                if (Nparm GE 6) then  pderiv(0,5) = 1
           endif

        if (Nparm GE 6) then return, factor * Lcf + parms(5) $
                        else return, factor * Lcf
     END

   else: BEGIN
        factor = parms(0,0)
        cntrd = parms(1,*)
        radius = parms(2,*)
        pscale = parms(3,*)
        power = parms(4,0)
        if (Nparm GE 6) then offset = parms(5,0) else offset=0
                return,0
           END

ENDCASE
end