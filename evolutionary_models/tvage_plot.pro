pro tvage_plot,output_dir=output_dir,model_name=model_name

cgcleanup

cgps_open,'/Users/gabeweible/Desktop/T_v_Age_'+model_name+'.ps',xsize=8,ysize=10,$
/portrait

;MWC 758

agec=3.5
tc=400
tcerrl=0
tcerrh=0

agecerr=2

plotsym,0,0.01,/fill

cgplot,[agec],[tc],err_xlow=[agecerr],err_xhigh=[agecerr],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,title=textoidl('T_{eff} vs. Age'),xtitle='Age (Myr)',$
ytitle=textoidl('T_{eff} (K)'),xrange=[.9,220],/xlog,yrange=[0,2700],charsize=2.5


alt_clr='grey'

;Marley+2008 cold start

folder='cold_start'
str='model_seq'
readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'.0010',n,t,l,r,ts,teff
;t is in yr, so let's go to Myr
t=t/1e6
;M is in Msun, so let's go to MJup
mark_offset=160
cgoplot,t,teff-mark_offset,color='dodger blue',linestyle=0
cgtext,140,220-mark_offset,textoidl('1 M_J'),color='dodger blue',size=1

folder='cold_start'
str='model_seq'
readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'.0020',n,t,l,r,ts,teff
;t is in yr, so let's go to Myr
t=t/1e6
;M is in Msun, so let's go to MJup
cgoplot,t,teff-mark_offset,color='dodger blue',linestyle=0
cgtext,140,300-mark_offset,textoidl('2 M_J'),color='dodger blue',size=1

folder='cold_start'
str='model_seq'
readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'.0040',n,t,l,r,ts,teff
;t is in yr, so let's go to Myr
t=t/1e6
;M is in Msun, so let's go to MJup
cgoplot,t,teff-mark_offset,color='dodger blue',linestyle=0
cgtext,140,390-mark_offset,textoidl('4 M_J'),color='dodger blue',size=1


folder='cold_start'
str='model_seq'
readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'.0080',n,t,l,r,ts,teff
;t is in yr, so let's go to Myr
t=t/1e6
;M is in Msun, so let's go to MJup
cgoplot,t,teff-mark_offset,color='dodger blue',linestyle=0,thick=8
cgtext,132,512-mark_offset,textoidl('8 M_J'),color='dodger blue',size=1

;folder='cold_start'
;str='model_seq'
;readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'.0080',n,t,l,r,ts,teff
;t is in yr, so let's go to Myr
;t=t/1e6
;M is in Msun, so let's go to MJup
;cgoplot,t,teff,color='dodger blue',linestyle=0
;cgtext,100,1600,textoidl('8 M_J'),color='dodger blue',size=0.7


folder='cold_start'
str='model_seq'
readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'.0100',n,t,l,r,ts,teff
;t is in yr, so let's go to Myr
t=t/1e6
;M is in Msun, so let's go to MJup
cgoplot,t,teff-mark_offset,color='dodger blue',linestyle=0
cgtext,130,660-mark_offset,textoidl('10 M_J'),color='dodger blue',size=1




;hot start

folder='ames_cond\ 3'
;folder='bt_settl'
;folder='bt_settl_agss2009_meta_05'
;folder='bt_dusty_agss2009_meta_00'

str='ames_cond'
;str='bt_settl'
;str='bt_settl_agss2009_meta_05'
;str='bt_dusty_agss2009_meta_00'

;evo tracks
readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'_m_0.0005.dat',t,m,teff
;t is in Gyr, so let's go to Myr
t=t*1000.
;M is in Msun, so let's go to MJup
m=m/0.00095
cgoplot,t,teff,color='crimson',linestyle=2
cgtext,1.2,650,textoidl('0.5 M_J'),color='crimson',size=1
readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'_m_0.001.dat',t,m,teff
;t is in Gyr, so let's go to Myr
t=t*1000.
;M is in Msun, so let's go to MJup
m=m/0.00095
cgoplot,t,teff,color='crimson',linestyle=2
cgtext,1.2,950,textoidl('1 M_J'),color='crimson',size=1

readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'_m_0.002.dat',t,m,teff
;t is in Gyr, so let's go to Myr
t=t*1000.
;M is in Msun, so let's go to MJup
m=m/0.00095
cgoplot,t,teff,color='crimson',linestyle=2
cgtext,1.2,1300,textoidl('2 M_J'),color='crimson',size=1

readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'_m_0.004.dat',t,m,teff
;t is in Gyr, so let's go to Myr
t=t*1000.
;M is in Msun, so let's go to MJup
m=m/0.00095
cgoplot,t,teff,color='crimson',linestyle=2
cgtext,1.2,1700,textoidl('4 M_J'),color='crimson',size=1

readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'_m_0.008.dat',t,m,teff
;t is in Gyr, so let's go to Myr
t=t*1000.
;M is in Msun, so let's go to MJup
m=m/0.00095
cgoplot,t,teff,color='crimson',linestyle=2,thick=8
cgtext,1.2,2000,textoidl('8 M_J'),color='crimson',size=1

readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'_m_0.01.dat',t,m,teff
;t is in Gyr, so let's go to Myr
t=t*1000.
;M is in Msun, so let's go to MJup
m=m/0.00095
cgoplot,t,teff,color='crimson',linestyle=2
cgtext,1.2,2220,textoidl('10 M_J'),color='crimson',size=1
readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'_m_0.012.dat',t,m,teff
;t is in Gyr, so let's go to Myr
t=t*1000.
;M is in Msun, so let's go to MJup
m=m/0.00095
cgoplot,t,teff,color='crimson',linestyle=2
cgtext,1.2,2300,textoidl('12 M_J'),color='crimson',size=1
readcol,'/Users/gabeweible/Downloads/'+folder+'/'+str+'_m_0.02.dat',t,m,teff
;t is in Gyr, so let's go to Myr
t=t*1000.
;M is in Msun, so let's go to MJup
m=m/0.00095
cgoplot,t,teff,color='crimson',linestyle=2
cgtext,1.2,2400,textoidl('20 M_J'),color='crimson',size=1

if 0 then begin
;CS prediction

plotsym,0,2,/fill

cgoplot,[agec],[775],color='black',psym=8
plotsym,0,1.25,/fill

cgoplot,[agec],[775],color='dodger blue',psym=8

;HS prediction
plotsym,0,2,/fill

cgoplot,[agec],[1900],color='black',psym=8
plotsym,0,1.25,/fill

cgoplot,[agec],[1900],color='crimson',psym=8
endif




;HR 8799b

agec=30
agecerrl=10
agecerrh=20

tc=870
tcerrl=70
tcerrh=30

clr=alt_clr

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr

cgtext,35,950,'HR8799bcde',charsize=1,color=alt_clr

;HR 8799c/d

agec=30
agecerrl=10
agecerrh=20

tc=1090
tcerrl=90
tcerrh=10

;clr='purple'

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr



;HR 8799e

agec=30
agecerrl=10
agecerrh=20

tc=1150
tcerrl=50
tcerrh=50

;clr='purple'

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr




;PDS 70b

agec=5.4
agecerrl=1
agecerrh=1

tc=1218
tcerrl=64
tcerrh=112

;clr='crimson'

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr


cgtext,6,1100,'PDS70bc',charsize=1,color=alt_clr


;PDS 70c

agec=5.4
agecerrl=1
agecerrh=1

tc=1030
tcerrl=289
tcerrh=216

;clr='crimson'

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr





;51 Eri b

agec=20
agecerrl=5
agecerrh=5

tc=700
tcerrl=100
tcerrh=50

clr=alt_clr

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr

cgtext,25,600,'51 Eri b',charsize=1,color=alt_clr



;AB Aur b

agec=3
agecerrl=2
agecerrh=2

tc=2200
tcerrl=0
tcerrh=0

clr=alt_clr

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr

cgtext,3,2250,'AB Aur b',charsize=1,color=alt_clr


;HIP75056Ab

agec=12
agecerrl=5
agecerrh=5

tc=2300
tcerrl=300
tcerrh=300

clr=alt_clr

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr

cgtext,15,2400,'HIP 75056Ab',charsize=1,color=alt_clr



;HIP65426b

agec=17
agecerrl=3
agecerrh=5

tc=1450
tcerrl=150
tcerrh=150

clr=alt_clr

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr

cgtext,6,1500,'HIP 65426b',charsize=1,color=alt_clr



;Beta Pic b
;https://arxiv.org/pdf/1703.00011.pdf
agec=24
agecerrl=3
agecerrh=3
;24 pm3, Bell et al. 2015
tc=1750
tcerrl=50
tcerrh=50

clr=alt_clr

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr

cgtext,12,1800,textoidl('\beta Pic b'),charsize=1,color=alt_clr



;HD 95086
;de Rosa 2016
agec=17
agecerrl=3
agecerrh=5
tc=1050
tcerrl=250
tcerrh=250

clr=alt_clr

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr

cgtext,8,880,textoidl('HD 95086b'),charsize=1,color=alt_clr



;AB Pic b
;bonnefoy 2010
agec=30
agecerrl=10
agecerrh=10
tc=2000
tcerrl=100
tcerrh=300

clr=alt_clr

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr

cgtext,50,2000,textoidl('AB Pic b'),charsize=1,color=alt_clr


;Kap And b
;Stone 2021
agec=42
agecerrl=6
agecerrh=4
tc=1800
tcerrl=300
tcerrh=300

clr=alt_clr

plotsym,0,1.5,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8

plotsym,0,1,/fill
cgoplot,[agec],[tc],err_xlow=[agecerrl],err_xhigh=[agecerrh],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8,color=clr

cgtext,50,1750,textoidl('\kappa And b'),charsize=1,color=alt_clr



;MWC 758 again

agec=3.5
tc=400
tcerrl=0
tcerrh=0

agecerr=2

plotsym,0,.01,/fill

cgoplot,[agec],[tc],err_xlow=[agecerr],err_xhigh=[agecerr],err_ylow=[tcerrl],$
err_yhigh=[tcerrh],psym=8



restore,'/Users/gabeweible/Desktop/'+output_dir+'/'+model_name+'_1.0RJup_avs_temps.sav'
;temps=reverse(temps)
;avs=reverse(avs)

;avs=reverse(avs)
;temps=reverse(temps)

temps=temps[where(avs mod 20 eq 0 and avs le 100)]
avs=avs[where(avs mod 20 eq 0 and avs le 100)]

avs=avs[uniq(temps)]
temps=temps[uniq(temps)]

for nn=0,n_elements(avs)-1 do begin
	;if avs[nn] mod 20 eq 0 then begin
	tcerrh=temps[nn]-400
	cgoplot,[agec],[tc],err_xlow=[agecerr],err_xhigh=[agecerr],err_ylow=[tcerrl],$
	err_yhigh=[tcerrh],color='black',err_thick=8
	if nn eq 0 then cgtext,agec+0.15,temps[nn]+20,textoidl('A_V='+sigfig(avs[nn],2)),$
	charsize=0.9
	if nn gt 0 then cgtext,agec+0.2,temps[nn]-30,textoidl('A_V='+sigfig(avs[nn],2)),$
	charsize=1
	;endif
endfor


restore,'/Users/gabeweible/Desktop/'+output_dir+'/'+model_name+'_2.0RJup_avs_temps.sav'
;temps=reverse(temps)
;avs=reverse(avs)

;avs=reverse(avs)
;temps=reverse(temps)


temps=temps[where(avs mod 20 eq 0)]
avs=avs[where(avs mod 20 eq 0)]

avs=avs[uniq(temps)]
temps=temps[uniq(temps)]

for nn=1,n_elements(avs)-2 do begin
	;if avs[nn] mod 20 eq 0 then begin
	tcerrh=temps[nn]-400
	cgoplot,[agec],[tc],err_xlow=[agecerr],err_xhigh=[agecerr],err_ylow=[tcerrl],$
	err_yhigh=[tcerrh],color='black',err_thick=8
	if nn le 2 then cgtext,agec-1.4,temps[nn]-30,textoidl('A_V='+sigfig(avs[nn],2)),$
	charsize=1
	if nn gt 2 then cgtext,agec-1.5,temps[nn]-30,textoidl('A_V='+sigfig(avs[nn],2)),$
	charsize=1
	;endif
endfor
	if nn gt 0 then cgtext,agec+0.2,330,textoidl('1 R_{Jup}'),charsize=1.1;,color='dodger blue'

	if nn gt 0 then cgtext,agec-1.3,330,textoidl('2 R_{Jup}'),charsize=1.1;,color='crimson'


cgtext,2.15,220,'MWC 758c',charsize=1.4,color='black'


cgps_close,/pdf





end