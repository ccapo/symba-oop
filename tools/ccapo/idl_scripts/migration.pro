pro migration, mass, data, migrate, res, fit, nofit=nofit
;
if(not keyword_set(nofit)) then nres = 2 else nres = 1
nfiles = n_elements(data) & nt = n_elements(data[0].t) & nsim = 2 & nsmooth = 2 & nm = n_elements(mass) & ntype = 3
migrate = replicate({a:dblarr(nres), a_err:dblarr(nres), e:dblarr(nres), e_err:dblarr(nres), i:dblarr(nres), i_err:dblarr(nres), m:dblarr(nres), m_err:dblarr(nres)}, nfiles/nsim)

; Loop over the different types of simulations for averaging purposes
for j = 0, nfiles/nsim - 1 do begin
B:
  nn = lonarr(nsim)
  if(j gt nfiles/nsim - 1) then goto, C

  for i = 0, nsim - 1 do begin
    dj = where(data[nsim*j + i].a gt 0.0d0, nj)
    nn[i] = nj
  endfor
  nd = max(nn,k)
  d = where(data[nsim*j + k].a gt 0.0d0)

  deltat = max(data[nsim*j].t[d]) - min(data[nsim*j].t[d]) ; Time span for integration
  slen = 2^nsmooth ; Smoothing length

  a = dblarr(nd, nsim) & as = dblarr(nd, nsim) & da = dblarr(nd, nsim) & das = dblarr(nd, nsim) & adot_rand = dblarr(nsim)
  ;e = dblarr(nd, nsim) & es = dblarr(nd, nsim) & de = dblarr(nd, nsim) & des = dblarr(nd, nsim) & edot_rand = dblarr(nsim)
  ;i = dblarr(nd, nsim) & is = dblarr(nd, nsim) & di = dblarr(nd, nsim) & dis = dblarr(nd, nsim) & idot_rand = dblarr(nsim)
  ;m = dblarr(nd, nsim) & ms = dblarr(nd, nsim) & dm = dblarr(nd, nsim) & dms = dblarr(nd, nsim) & mdot_rand = dblarr(nsim)
  for i = 0, nsim - 1 do begin
    ; Store data for conveinence
    a[*,i] = data[nsim*j + i].a[d]
    ;e[*,i] = data[nsim*j + i].e[d]
    ;i[*,i] = data[nsim*j + i].i[d]
    ;m[*,i] = data[nsim*j + i].mass[d]

    ; Spacing between adjacent data points
    da[*,i] = shift(a[*,i],-1) - a[*,i]
    ;de[*,i] = shift(e[*,i],-1) - e[*,i]
    ;di[*,i] = shift(i[*,i],-1) - i[*,i]
    ;dm[*,i] = shift(m[*,i],-1) - m[*,i]

    ; Mean migration rate assuming a 1-D random walk
    adot_rand[i] = stddev(da[0:nd - 2,i])*sqrt(nd)/deltat
    ;edot_rand[i] = stddev(de[0:nd - 2,i])*sqrt(nd)/deltat
    ;idot_rand[i] = stddev(di[0:nd - 2,i])*sqrt(nd)/deltat
    ;mdot_rand[i] = stddev(dm[0:nd - 2,i])*sqrt(nd)/deltat

    ; Smooth data to improve probability of finding a large continous monotonic sequence
    as[*,i] = smooth(a[*,i],slen,/edge)
    ;es[*,i] = smooth(e[*,0],slen,/edge)
    ;is[*,i] = smooth(i[*,0],slen,/edge)
    ;ms[*,i] = smooth(m[*,0],slen,/edge)

    ; Spacing between adjacent smoothed data points
    das[*,i] = shift(as[*,i],-1) - as[*,i]
    ;des[*,i] = shift(es[*,i],-1) - es[*,i]
    ;dis[*,i] = shift(is[*,i],-1) - is[*,i]
    ;dms[*,i] = shift(ms[*,i],-1) - ms[*,i]
  endfor

  ; Find the slope of the linear migration regime
  if(not keyword_set(nofit)) then begin
    adot_lin = fit_data(data[nsim*j].t[d[0:nd - 2]], das[0:nd - 2,*], a[0:nd - 2,*])
    ;edot_lin = fit_data(data[nsim*j].t[d[0:nd - 2]], des[0:nd - 2,*], e[0:nd - 2,*])
    ;idot_lin = fit_data(data[nsim*j].t[d[0:nd - 2]], dis[0:nd - 2,*], i[0:nd - 2,*])
    ;mdot_lin = fit_data(data[nsim*j].t[d[0:nd - 2]], dms[0:nd - 2,*], m[0:nd - 2,*])

    ; Mean migration by fitting to the linear migration regime
    migrate[j].a[0] = mean(adot_lin)
    migrate[j].a_err[0] = stddev(adot_lin)

    ;migrate[j].e[0] = mean(edot_lin)
    ;migrate[j].e_err[0] = stddev(edot_lin)

    ;migrate[j].i[0] = mean(idot_lin)
    ;migrate[j].i_err[0] = stddev(idot_lin)

    ;migrate[j].m[0] = mean(mdot_lin)
    ;migrate[j].m_err[0] = stddev(mdot_lin)
  endif

  ; Mean migration assuming a 1-D random walk
  migrate[j].a[nres - 1] = mean(adot_rand)
  migrate[j].a_err[nres - 1] = stddev(adot_rand)

  ;migrate[j].e[nres - 1] = mean(edot_rand)
  ;migrate[j].e_err[nres - 1] = stddev(edot_rand)

  ;migrate[j].i[nres - 1] = mean(idot_rand)
  ;migrate[j].i_err[nres - 1] = stddev(idot_rand)

  ;migrate[j].m[nres - 1] = mean(mdot_rand)
  ;migrate[j].m_err[nres - 1] = stddev(mdot_rand)
endfor
C:

res = replicate({a:dblarr(2, ntype), e:dblarr(2, ntype), i:dblarr(2, ntype), m:dblarr(2, ntype)}, nres)
fit = replicate({a:dblarr(nm, ntype), e:dblarr(nm, ntype), i:dblarr(nm, ntype), m:dblarr(nm, ntype)}, nres)
for i = 0, nres - 1 do begin

  for j = 0, ntype - 1 do begin

    ; index range
    kmin = j*nm
    kmax = (j + 1)*nm - 1

    ; fits of adot vs mass
    d = where(migrate[kmin:kmax].a_err[i] ne 0.0d0, dcount)
    if(dcount eq n_elements(mass)) then begin
      res[i].a[*,j] = linfit(mass, migrate[kmin:kmax].a[i], measure_errors = migrate[kmin:kmax].a_err[i], yfit = tmp, /double)
      fit[i].a[*,j] = tmp
    endif else begin
      res[i].a[*,j] = linfit(mass, migrate[kmin:kmax].a[i], yfit = tmp, /double)
      fit[i].a[*,j] = tmp
    endelse

    ; fits of edot vs mass
    ;d = where(migrate[kmin:kmax].e_err[i] ne 0.0d0, dcount)
    ;if(dcount eq n_elements(mass)) then begin
    ;  res[i].e[*,j] = linfit(mass, migrate[kmin:kmax].e[i], measure_errors = migrate[kmin:kmax].e_err[i], yfit = tmp, /double)
    ;  fit[i].e[*,j] = tmp
    ;endif else begin
    ;  res[i].e[*,j] = linfit(mass, migrate[kmin:kmax].e[i], yfit = tmp, /double)
    ;  fit[i].e[*,j] = tmp
    ;endelse

    ; fits of idot vs mass
    ;d = where(migrate[kmin:kmax].i_err[i] ne 0.0d0, dcount)
    ;if(dcount eq n_elements(mass)) then begin
    ;  res[i].i[*,j] = linfit(mass, migrate[kmin:kmax].i[i], measure_errors = migrate[kmin:kmax].i_err[i], yfit = tmp, /double)
    ;  fit[i].i[*,j] = tmp
    ;endif else begin
    ;  res[i].i[*,j] = linfit(mass, migrate[kmin:kmax].i[i], yfit = tmp, /double)
    ;  fit[i].i[*,j] = tmp
    ;endelse

    ; fits of mdot vs mass
    ;d = where(migrate[kmin:kmax].m_err[i] ne 0.0d0, dcount)
    ;if(dcount eq n_elements(mass)) then begin
    ;  res[i].m[*,j] = linfit(mass, migrate[kmin:kmax].m[i], measure_errors = migrate[kmin:kmax].m_err[i], yfit = tmp, /double)
    ;  fit[i].m[*,j] = tmp
    ;endif else begin
    ;  res[i].m[*,j] = linfit(mass, migrate[kmin:kmax].m[i], yfit = tmp, /double)
    ;  fit[i].m[*,j] = tmp
    ;endelse

  endfor

endfor

end

function fit_data, time, ddata, data
;
nt=n_elements(ddata[*,0]) & nsim = n_elements(ddata[0,*])
res = dblarr(2, nsim)
for j = 0, nsim - 1 do begin
  ; Search for longest continuous block of either monotonically increasing or decreasing data
  pos = where(ddata[*,j] ge 0.0d0,pcount) & neg = where(ddata[*,j] le 0.0d0,ncount)
  if(pcount ge 2) then begin
    dpos = shift(pos,-1)-pos & dpos = dpos[0:n_elements(dpos)-2]
    idpos = where(dpos ne 1,pcount2)
    if(pcount2 ne 0) then begin
      didpos = idpos-shift(idpos,1) & didpos[0] = idpos[0]
      maxpos = max(didpos,imaxpos)
    endif else begin
      idpos = [0, nt - 1] & maxpos = nt - 1 & imaxpos = 1
    endelse
  endif else maxpos = 0
  if(ncount ge 2) then begin
    dneg = shift(neg,-1)-neg & dneg = dneg[0:n_elements(dneg)-2]
    idneg = where(dneg ne 1,ncount2)
    if(ncount2 ne 0) then begin
      didneg = idneg-shift(idneg,1) & didneg[0] = idneg[0]
      maxneg = max(didneg,imaxneg)
    endif else begin
      idneg = [0, nt - 1] & maxneg = nt - 1 & imaxneg = 1
    endelse
  endif else maxneg = 0
  if(maxpos ne maxneg) then begin
    maxid = max([maxpos, maxneg],imaxid)
    if(imaxid eq 0) then begin
      if(imaxpos eq 0) then idmin = 0 else idmin = idpos[imaxpos - 1]
      idmax = idpos[imaxpos]
    endif else begin
      if(imaxneg eq 0) then idmin = 0 else idmin = idneg[imaxneg - 1]
      idmax = idneg[imaxneg]
    endelse
  endif else begin
    idmin = 0 & idmax = nt - 1
  endelse

  ; Check if selected portion is alright
A:
  ans = ''
  id = idmin + indgen(idmax - idmin)
  plot, time, data[*,j] & oplot, time[id], data[id,j], color = 254
  print, 'Is the desired portion highlighted?'
  read, ans
  if(ans eq 'n') then begin
    print, 'Enter values for the index range (i.e. idmin,idmax) '
    read, idmin, idmax
    goto, A
  endif

  ; Now perform the fit
  ans = ''
  res[*,j] = linfit(time[id], data[id,j], yfit = yfit, /double)
  plot, time[id], data[id,j] & oplot, time[id], yfit, line = 2
  print, 'Does the fit look alright?'
  read, ans
  if(ans eq 'n') then goto, A
endfor

return, reform(res[1,*])
end
