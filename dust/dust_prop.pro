pro dust_prop,type

; Possible valuse for type (all STRING type):
; - g
; - albedo
; - kappa

base='astrosil_'
gsize=['0.1','0.15','0.2','0.3','0.4','0.5','0.75','1','1.5','2','3','4','5','7.5','10','15','20','30','40','50',$
      '75','100','150','200','300','400','500','750','1000','1500','2000','3000']
ndir=n_elements(gsize)

wvl=readfits(strcompress(base+gsize(0)+'mic/lambda.fits.gz',/remove_all),/silent)
nwvl=n_elements(wvl)

props=fltarr(nwvl,ndir)

for i=0,ndir-1 do begin
   prop=readfits(strcompress(base+gsize(i)+'mic/'+type+'.fits.gz',/remove_all),/silent)
   props(*,i)=prop
endfor

xmin=min(wvl)/1.2
xmax=max(wvl)*1.2
case type of
   'g': begin
      ymax=1.
      ymin=0.
      plot,[0.],[1e-3],/nodata,xrange=[xmin,xmax],yrange=[ymin,ymax],/xstyle,/ystyle,/xlog,$
           xtitle='!7k (l!6m)',ytitle='!6g'
   end
   'albedo': begin
      ymax=1.
      ymin=0.
      plot,[0.],[1e-3],/nodata,xrange=[xmin,xmax],yrange=[ymin,ymax],/xstyle,/ystyle,/xlog,$
           xtitle='!7k (l!6m)',ytitle='!6Albedo'
   end
   'kappa': begin
      ymax=max(props)*1.2
      ymin=min(props)/1.2
      plot,[0.],[1e-3],/nodata,xrange=[xmin,xmax],yrange=[ymin,ymax],/xstyle,/ystyle,/xlog,$
           xtitle='!7k (l!6m)',ytitle='!7j!6 (cm!u2!n/g)',/ylog
   end
endcase

cols=findgen(ndir)*(255/ndir)

loadct,13
for i=0,ndir-1 do begin
   oplot,wvl,props(*,i),color=cols(i)
endfor
loadct,0

end
