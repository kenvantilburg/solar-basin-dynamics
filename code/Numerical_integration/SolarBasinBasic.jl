


using LinearAlgebra, StaticArrays

const year = 2pi
const rsun = 0.00463673;

# solar potential including potential inside the sun (analytical approximation).  rebound only deals in 
# 1/r potentials, and so we modify here.
sun_phi_approx(r) = 0.0019858674079720195 + 0.5010980317741527*((10.284885633196014 + 2.0/r)*exp(-10.284885633196014*r) - 2.0/r)

# further processing of the solar potential
function sun_phi(r)
	rr = r/rsun
	if rr < 1
		return sun_phi_approx(rr)/rsun
	end
	return -1/r
end

function rmin_kepler(x,v)
	e = 0.5*dot(v,v) - 1/norm(x)
	l = cross(x,v)
	l2 = dot(l,l)
	(-1 + sqrt(1 + 2*e*l2))/(2*e)
end

function r1xv(x,v)
	l = cross(x,v)
	ev = cross(v,l) - x/norm(x)
	ecc = norm(ev)
	e = 0.5*dot(v,v) - 1/norm(x)
	a = -0.5/e
	if a < 0. || ecc > 1.
        return nothing
    end
	bv = cross(l,ev)
	bv /= norm(bv)
	bv *= a*sqrt(1-ecc^2)
	av = a*ev/ecc
	#
	cn = (1-1/a)/ecc
	if abs(cn) > 1
		return nothing
	end
	sn = sqrt(1 - cn^2)
	dndt = 1/(a^1.5 * (1 - ecc*cn))
	acb = cross(av,bv)
	sgn = dot(l,acb) > 0. ? 1. : -1.
	return (((-ecc + cn)*av  + sn*bv, sgn*dndt*(-sn*av + cn*bv)), 
		   ((-ecc + cn)*av  - sn*bv, -sgn*dndt*(-sn*av - cn*bv)))
end
