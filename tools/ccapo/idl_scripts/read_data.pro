pro read_data, files, template, data, mass
;
nfiles = n_elements(files.name) & nt = 101 & nsim = 2
data = replicate({t:dblarr(nt), a:dblarr(nt), e:dblarr(nt), i:dblarr(nt), mass:dblarr(nt)}, nfiles)
mass = dblarr(nfiles)

; Import data files
for i = 0, nfiles - 1 do begin
  tmp = read_ascii(files.name[i]+'.out',template=template)
  ; Apparently all the steps below are necessary
  data[i].t = tmp.t
  data[i].a = tmp.a
  data[i].e = tmp.e
  data[i].i = tmp.i
  data[i].mass = tmp.mass
  mass[i] = data[i].mass[0]
endfor
mass = mass[sort(mass)]
mass = mass[uniq(mass)]

end
