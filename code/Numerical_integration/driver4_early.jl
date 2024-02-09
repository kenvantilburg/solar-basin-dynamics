# runs simulations for specific "early" particle indices

using OrdinaryDiffEq, LinearAlgebra, Random, StaticArrays, Printf

const npl = 4
const nplanets  = 5
const mp = @SVector [2.4478382877847715e-06, 3.040432648022642e-06, 0.0009547919152112404, 0.0002858856727222417, 0.0]

const x0p = [[-0.722774,-0.0196495,0.041168] [0.349701,0.92679,1.88066e-05] [2.84897,-4.2281,-0.0461854] [5.35027,-8.42927,-0.0663824] [0.,0.,0.]]

const v0p = [[-0.0352425,1.18056,0.0182419] [0.950192,-0.35493,1.54367e-06] [-0.358477,-0.265794,0.0091223] [-0.255767,-0.172943,0.0131922] [0.,0.,0.]]


const year = 2pi
const rsun = 0.00463673;
sun_mf_approx(r) = 0.5010980317741527*(2 + (-2 - 10.284885633196014*r*(2 + 10.284885633196014*r)) * exp(-10.284885633196014*r))

#################################################

function qdot!(dq,p,q,params,t)
	dq .= p
end

function v3(q,i)
	return @inbounds SVector(q[1,i],q[2,i],q[3,i])
end

function sv3(q,i,v)
	@inbounds begin
		q[1,i] += v[1]
		q[2,i] += v[2]
		q[3,i] += v[3]
	end
end

# Using MutableArray?
# Doing something more intelligent with test particles?
function pdot!(dp,p,q,params,t)
	# Is this necessary?
	fill!(dp,0.)
	x3 = @SVector [norm(v3(q,i))^3 for i in 1:nplanets]
	for i in 1:nplanets
		xi = v3(q,i)
		@inbounds xi3 = x3[i]
		for j in i+1:nplanets
			xj = v3(q,j)
			# Recalc ...
			@inbounds xj3 = x3[j]
			dx = xi .- xj
			dx3 = norm(dx)^3
			@inbounds sv3(dp, i, -mp[j]*(dx/dx3 + xj/xj3))
			@inbounds sv3(dp, j, mp[i]*(dx/dx3 - xi/xi3))
		end
	end
	# Recalc ...
	for i in 1:nplanets
		xi = v3(q,i)
		xi1 = norm(xi)
		#xi3 = xi1^3
		@inbounds xi3 = x3[i]
		@inbounds pre = 1 + mp[i]
		if xi1 < rsun
			pre *= sun_mf_approx(xi1/rsun)
		end
		sv3(dp,i, -pre*xi/xi3)
	end
end

#################################################

# TODO - checking effect of perturbations on this ...
function inject_outside_sun(rinj,arange)
	mf = 1.
	if rinj < rsun
		mf = sun_mf_approx(rinj/rsun)
	end
	vmin = sqrt(2*mf/rinj-1/arange[1])
	vmax = sqrt(2*mf/rinj-1/arange[2])
	vemit = (vmin^3 + (vmax^3-vmin^3)*rand())^(1/3)
	ul = rand()
	l0 = rinj*vemit*sqrt(1 - ul^2)
	e0 = 0.5*vemit^2 - mf/rinj
	a = -0.5/e0
	ecc = sqrt(1 + 2*e0*l0^2)
	b = a * sqrt(1-ecc^2)
	#
	av = randn(3)
	av = (a/norm(av)).*av
	# Then, random perpendicular vector
	zv = randn(3)
	bv = cross(zv,av)
	bv = (b/norm(bv)).*bv
	#
	eta = pi - (1-ecc)*(rand()-0.5)
	xv = (ecc-cos(eta)).*av .+ sin(eta).*bv
	vv = sin(eta).*av .+ cos(eta).*bv
	vr = sqrt(2*(e0 + 1/norm(xv)))
	vv = (vr/norm(vv)).*vv
	#
	(SVector{3}(xv),SVector{3}(vv))
end

function inject_sun(rinj,arange)
	mf = 1.
	if rinj < rsun
		mf = sun_mf_approx(rinj/rsun)
	end
	vmin = sqrt(2*mf/rinj-1/arange[1])
	vmax = sqrt(2*mf/rinj-1/arange[2])
	vemit = (vmin^3 + (vmax^3-vmin^3)*rand())^(1/3)
	# Then, choosing thing ...
	x = @SVector randn(3)
	x *= rinj/norm(x)
	v = @SVector randn(3)
	v *= vemit/norm(v)
	if dot(x,v) < 0.
		v = -v
	end
	#return (x,v)
	x0p[:,5] .= x
	v0p[:,5] .= v
end

function grav_pot(x)
	p = 0.
	for i in 1:4
		p -= mp[i]/norm(x - x0p[:,i])
	end
	p -= 1/norm(x)
	p
end

# Should be slight correction here from non-inertial
# frame ... fixing ...

# Assume Earth position in x0p ...
function inject_around_earth(de=0.002)
	xe = x0p[:,2]
	dx = @SVector randn(3)
	dx *= de/norm(dx)
	x = xe + dx
	# Then, choosing velocity ...
	gp = grav_pot(x)
	vmax = sqrt(-2*gp)
	ul = rand()
	v = @SVector randn(3)
	v *= ul^(1/3) * vmax/norm(v)
	# Then, setting ...
	x0p[:,5] .= x
	v0p[:,5] .= v
end

#################################################
#################################################


function printout(io,t0,tstep)
	@printf(io,"%a,%a",t0,tstep)
	for i in 1:nplanets
		@printf(io,",%a,%a,%a,%a,%a,%a",x0p[1,i],x0p[2,i],x0p[3,i],
				v0p[1,i],v0p[2,i],v0p[3,i])
	end
	@printf(io,"\n")
end

# Assumes that x0p, v0p have been set
# Appends to given file
function integ_file(outfile,t0,tgap,tmax,tstep,print_first)

	tspan = (0,1.5*tgap)

	io = open(outfile,"a")

	if print_first
		printout(io,t0,tstep)
	end

	if t0 > tmax || norm(v3(x0p,5)) > 30.
		close(io)
		return
	end

	i = 0

	while true
		prob = DynamicalODEProblem(pdot!,qdot!,v0p,x0p,tspan)
		integrator = init(prob,DPRKN12(),reltol=1e-14,abstol=1e-14,dt=tstep,maxiters=Inf)

		while true
			step!(integrator)
			err = check_error(integrator)
			if err != :Success
				println(err)
				close(io)
				return
			end

			if integrator.t > tgap
				t0::Float64 = t0 + integrator.t
				tstep::Float64 = integrator.t - integrator.tprev
				x0p .= integrator.u[2,:,:]
				v0p .= integrator.u[1,:,:]
				break
			end
		end

		printout(io,t0,tstep)

		if t0 > tmax || norm(v3(x0p,5)) > 30.
			break
		end

		if i % 50 == 0
			flush(io)
		end
		i += 1
	end

	flush(io)
	close(io)
end

function load_planets(i,earth)
	a = zeros(2 + npl*6)
	
	file = (earth ? 
			"starts_earth/start_$i.csv" :
			"starts_sun/start_$i.csv")

	open(file, "r") do f
		l = readline(f)
		nums = split(l,",")
		a = parse.(Float64,nums)
	end

	t0 = a[1]
	tstep = a[2]
	for i in 1:npl
		i0 = 2 + 6*(i-1)
		x0p[:,i] .= a[i0+1:i0+3]
		v0p[:,i] .= a[i0+4:i0+6]
	end

	return t0
end

function extend_file(outfile,tmax,tgap)
	a = zeros(2 + nplanets*6)

	open(outfile) do f
		seekend(f)
		pos1 = position(f)
		if pos1 >= 1000
			skip(f,-1000)
		else
			seek(f,0)
		end
		l = readline(f)
		if position(f) != pos1
			l = readline(f)
		end
		nums = split(l,",")
		a .= parse.(Float64,nums)
	end

	t0 = a[1]
	tstep = a[2]
	for i in 1:nplanets
		i0 = 2 + 6*(i-1)
		x0p[:,i] .= a[i0+1:i0+3]
		v0p[:,i] .= a[i0+4:i0+6]
	end

	integ_file(outfile,t0,tgap,tmax,tstep,false)
end

early_inds = [1,4,5,10,11,12,14,21,28,30,31,37,39,50,52,53,56,57,58,59,62,63,68,69,71,74,81,83,85,87,88,90,94,95,98,101,119,123,126,129,133,137,139,144,146,148,157,158,166,169,170,171,172,174,176,180,183,184,186,196,198,201,203,204,212,214,218,221,225,231,233,240,241,242,243,244,245,247,251,254,257,272,273,274,275,277,278,281,283,284,285,287,288,289,291,292,297,304,305,307,311,316,322,323,326,334,338,339,340,341,345,352,356,358,360,368,369,370,371,372,376,379,383,385,386,389,393,394,396,398,399,400,401,402,403,411,417,420,422,423,426,429,431,435,436,440,441,443,445,448,451,453,454,455,456,457,459,462,464,466,468,469,470,475,476,477,478,479,482,485,486,487,489,491,493,496,500,502,503,504,511,512,514,516,520,521,523,531,532,534,536,537,538,539,540,546,547,548,549,554,557,559,565,566,569,572,575,577,585,588,590,596,601,603,604,606,611,612,617,624,627,630,637,639,640,643,647,648,650,652,653,658,664,665,669,671,680,681,689,690,694,695,696,699,700,702,705,707,709,713,716,720,721,725,728,729,730,731,733,736,745,748,749,751,752,754,755,756,761,762,763,769,775,777,780,784,788,790,791,792,800,808,810,811,814,815,819,821,823,828,829,831,834,835,837,838,845,846,851,853,855,856,863,866,870,871,872,873,874,875,877,878,879,880,883,884,885,886,887,888,894,895,899,901,902,903,904,907,909,910,912,913,915,917,919,926,927,928,932,933,935,936,938,941,947,948,950,951,953,958,959,961,963,964,966,971,972,976,977,981,982,984,986,987,988,991,997,999]

# Loop function ...
function generate_files(i,outbase,seed,tmax,tgap,earth)
	while i < 1000
		si = seed + i*1000
		Random.seed!(si)
		# Needing to keep track of where we are ...
		open("$(outbase)_$seed/its","a") do f
			println(f,i)
		end
		#
		#planets_i = rand(1:1000)
		planets_i = rand(early_inds)
		outfile = "$(outbase)_$seed/output_$si.csv"
		# Then, generate initial data ...
		if earth
			t0 = load_planets(planets_i,true)
			inject_around_earth()
			integ_file(outfile,t0,tgap,tmax,2e-4,true)
		else
			t0 = load_planets(planets_i,false)
			inject_sun(rsun,(0.4,4.0))
			integ_file(outfile,t0,tgap,tmax,2e-5,true)
		end
		i += 1
	end
end

function run_script()
	if length(ARGS) == 0
		return
	end

	seed = parse(Int,ARGS[2])
	ty = parse(Float64,ARGS[3])
	tyg = parse(Float64,ARGS[4])
	outbase = ARGS[5]

	tmax = ty*year
	tgap = tyg*year

	earth_inj = false
	if ARGS[1] == "EARTH"
		earth_inj = true
	elseif ARGS[1] != "SUN"
		# Something failed here ...
		return
	end

    # See if there's existing stuff ...
    ilast = 0
	its_file = "$(outbase)_$seed/its"
	if isfile(its_file)
		open(its_file,"r") do f
			for line in readlines(f)
				ilast = parse(Int,line)
			end
		end
		si = seed + ilast*1000
		outfile = "$(outbase)_$seed/output_$si.csv"
		if isfile(outfile)
			extend_file(outfile,tmax,tgap)
			ilast += 1
		end
	end
	# Start generating new files, after old one has been extended
	generate_files(ilast,outbase,seed,tmax,tgap,earth_inj)
end

run_script()
