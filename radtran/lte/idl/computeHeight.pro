; Compute the height axis from the opacity at 500 nm and the optical depth at the same wavelength
function computeHeight, tau500, opacity500
	n = n_elements(tau500)
	z = fltarr(n)
	
	z[0] = 0.d0
	for k = 1, n-1 do begin
		avgOpacity = (opacity500[k] + opacity500[k-1]) / 2.d0
		z[k] = z[k-1] - (tau500[k] - tau500[k-1]) / avgOpacity
	endfor
	
; Put z=0 at the point where tau500=1
	res = min(abs(tau500 - 1.d0), loc)
	z = z - z[loc]
	
	return, z
end