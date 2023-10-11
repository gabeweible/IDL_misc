pro tidal_decay_full

output=1


M1 = double(0.6);E30
M2 = double(0.6);E30
Msun=double(2.0E30)
R1 = double(0.5*695700./(1.5E8));
ain = double(0.44);0.12);
ein = double(0.5);
;Period = double(9.93);days
G=double(6.67E-11)
Einitial = - double(M1*M2/(2.*ain));scaled
newe = Einitial;
a = ain;
ecc = ein;
rperi=ain*(1.d0-ein)
time=0.d0

omega1=double(2.d0*!PI*1.7*695700./150./86400.)
omegaAB=double(9.93)
period=omegaAB

nn=0.

Ftide=20.d0	;check this ater
k1=double(0.13 - (0.11) *(alog(14.)/alog(27.)))
tau=double(0.35 - (0.35-0.05)*(alog(14.)/alog(27.)))	;in years
q2=double(M2/M1)
tau=double(tau*365.25) ;keeping tau in days


timestep=1.

;dE = double(0.3*2.16*(M1 + M2)/M1  * ((M2)^2.)/R1 *((a *(1. - ecc))/R1)^(-9.));


for ii=0., 3.E7,1. do begin	;currently in orbits, need to change to timesteps at least on the order of a day
						;UPDATE: now is in days
dt=timestep

if double(time)/365.25 lt 16.0E6  and a gt 0 and ecc gt 0 then begin

;newe = double(1.0*newe) - double(dE);




;if newe eq 32. then newe=33.
;a = ain*Einitial/newe;

dA= dt * ( -6.d0 *  Ftide * (k1/tau) * q2 * (1.d0+q2)*( (R1/a)^(8.d0) ) * a/((1.d0-ecc*ecc)^(15.d0/2.d0))) * $
	(1.d0 + (31.d0/2.d0)*ecc*ecc + (255.d0/8.d0)*ecc^4.d0 + (185.d0/16.d0)*ecc^6.d0 + (25.d0/64.d0)*ecc^8.d0 $
	-(1.d0-ecc*ecc)^(3.d0/2.d0) * (1.d0+ (15.d0/2.d0)*ecc^2.d0 + (45.d0/8.d0)*ecc^4.d0 + (5.d0/16.d0)*ecc^6.d0 )* (omega1/omegaAB)) ;(1.d0+6.d0*ecc^2.d0+(3.d0/8.d0)*ecc^4.d0+(233.d0/8.d0)*ecc^6.d0)) 
	
dE= dt * (-27.d0)*Ftide * (k1/tau) * q2 * (1.d0+q2)*( (R1/a)^(8.d0) ) * ecc/((1.d0-ecc*ecc)^(13.d0/2.d0)) * $
	(1.d0+(15.d0/4.d0)*ecc^2.d0+(15.d0/8.d0)*ecc^4.d0+(5.d0/64.d0)*ecc^6.d0 -(11.d0/8.d0)*((1.d0-ecc*ecc)^(3.d0/2.d0)) * (  1.d0+(3.d0/2.d0)*ecc^2.d0+(1.d0/8.d0)*ecc^4.d0  ) *(omega1/omegaAB) );(1.d0+6.d0*ecc^2.d0+(3.d0/8.d0)*ecc^4.d0+(233.d0/8.d0)*ecc^6.d0) )


sync=(1.d0+6.d0*ecc^2.d0+(3.d0/8.d0)*ecc^4.d0+(233.d0/8.d0)*ecc^6.d0)

a=a+dA

ecc=ecc+dE

newe=-double(M1*M2/(2.*a))

;ecc=1.-(rperi/a)
Period = Sqrt((a^3.)/(M1 + M2))*365.25;
omegaAB=Period
time=time+period

ecc=max([0.,ecc])
a=max([0.,a])

if time/365.25 ge nn*1E4 then begin 

if output then begin
print, '-----Finished tidal decay for timestep # ',ii

print, 'Initial E = ', Einitial

print, 'dA = ',dA

print, 'd = ', dE

print, 'Sync = ', sync

print, 'Period = ',period


print, 'Energy new = ', newe

print
;print, 'dE = ', dE
print
print, 'a = ',a
print
print, 'e = ',ecc
print
print, 'E = ',newe/Einitial*100.,'% initial value'
print
print, 'Time thus far = ',time/365.25,' years.'
endif

if nn eq 0 then begin
	times=time
	as=a
	eccs=ecc
	energies=newe
	syncs=sync
endif else begin
	times=[times,time]
	as=[as,a]
	eccs=[eccs,ecc]
	energies=[energies,newe]
	syncs=[syncs,sync]
endelse


nn=nn+1


endif

endif else begin


	;plot,times/365.5,eccs
	;oplot,times/365.5,as


print, '-----Finished tidal decay for timestep # ',ii

print, 'Initial E = ', Einitial

print, 'E new = ', newe

print, 'dA = ', dA

print, 'dE = ', dE
print
print, 'a = ',a
print
print, 'e = ',ecc
print
print, 'E = ',newe/Einitial*100.,'% initial value'
print
print, 'Time thus far = ',time/365.25,' years.'



;do plotting here

mydevice = !D.NAME
; Set plotting to PostScript:
SET_PLOT, 'PS'
; Use DEVICE to set some PostScript device options:
DEVICE, FILENAME='~/Desktop/AB.ps', /COLOR,/INCHES,/PORTRAIT;, YSIZE=8.
; Makes a new plot for every spectrum

!p.charsize=1.25
!p.thick=2


!p.multi=[0,3,2]


for kk=0,2 do begin

if kk eq 0 then yaxis=as
if kk eq 1 then yaxis=eccs
if kk eq 2 then yaxis=energies



if kk eq 0 then y_title='Semi-major axis (au)'
if kk eq 1 then y_title='Eccentricity'
if kk eq 2 then y_title=textoidl('Energy / | Initial Energy |')

x_title='Time (Myr)'


xaxis=times/365.25/1000000.

plot,xaxis,yaxis,xtitle=x_title,ytitle=y_title,xthick=3,ythick=3,charthick=3

tvlct,0,0,255,3
!p.color=3
!p.thick=5

oplot,xaxis,yaxis

tvlct,0,0,0,3

!p.color=0


endfor






DEVICE, /CLOSE
; Return plotting to the original device:
SET_PLOT, mydevice




return



endelse

endfor

beep





end