# driver4 but just runs a select few files

using OrdinaryDiffEq, LinearAlgebra, Random, StaticArrays, Printf

const npl = 4
const nplanets  = 5
const mp = @SVector [2.4478382877847715e-06, 3.040432648022642e-06, 0.0009547919152112404, 0.0002858856727222417, 0.0]

const x0p = [[-0.722774,-0.0196495,0.041168] [0.349701,0.92679,1.88066e-05] [2.84897,-4.2281,-0.0461854] [5.35027,-8.42927,-0.0663824] [0.,0.,0.]]

const v0p = [[-0.0352425,1.18056,0.0182419] [0.950192,-0.35493,1.54367e-06] [-0.358477,-0.265794,0.0091223] [-0.255767,-0.172943,0.0131922] [0.,0.,0.]]


const year = 2pi
const rsun = 0.00463673;
sun_mf_approx(r) = 0.5010980317741527*(2 + (-2 - 10.284885633196014*r*(2 + 10.284885633196014*r)) * exp(-10.284885633196014*r))

const priority_list = [88, 187, 43, 57, 29, 220, 138, 189, 72, 89, 234, 195, 24, 155, 147, 16, 18, 40, 246, 9, 77, 19, 175, 101, 20, 247, 255, 32, 47, 215, 134, 46, 205, 231, 196, 70, 180, 125, 60, 115, 199, 45, 105, 172, 64, 59, 2, 208, 217, 240, 78, 54, 28, 71, 252, 228, 143, 36, 148, 58, 17, 106, 197, 182, 190, 25, 139, 66, 49, 75, 211, 118, 100, 239, 119, 114, 87, 219, 224, 31, 193, 203, 243, 35, 97, 112, 109, 162, 166, 102, 128, 52, 113, 39, 186, 80, 202, 42, 178, 226, 131, 225, 26, 176, 8, 50, 223, 170, 91, 62, 63, 216, 111, 213, 233, 210, 221, 192, 167, 254, 145, 163, 21, 198, 83, 65, 34, 158, 185, 23, 110, 153, 13, 152, 30, 14, 73, 212, 67, 4, 149, 165, 242, 207, 103, 6, 256, 230, 194, 85, 142, 15, 74, 96, 129, 232, 251, 10, 95, 204, 5, 181, 7, 137, 81, 53, 68, 133, 206, 156, 61, 200, 241, 244, 86, 120, 171, 92, 222, 161, 136, 169, 38, 218, 229, 214, 55, 3, 51, 33, 124, 183, 168, 188, 238, 123, 132, 107, 93, 160, 179, 250, 154, 201, 253, 116, 127, 135, 94, 69, 104, 184, 76, 11, 173, 48, 98, 84, 144, 235, 117, 141, 108, 130, 41, 90, 126, 164, 157, 1, 122, 27, 151, 12, 227, 159, 249, 22, 79, 150, 209, 140, 146, 56, 236, 44, 237, 248, 245, 121, 177, 37, 174, 99, 82, 191]
const plist = priority_list[1:32]

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
		planets_i = rand(1:1000)
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

	# Check if we're in the priority list
	if !in(seed,plist)
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
end

run_script()

