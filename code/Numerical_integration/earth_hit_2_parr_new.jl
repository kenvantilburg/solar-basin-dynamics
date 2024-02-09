

#using Images
using JLD2
using Distributed

@everywhere using SharedArrays

include("OutputFiles2.jl")
@everywhere include("SolarBasinBasic.jl")

@everywhere using OrdinaryDiffEq

@everywhere const nplanets  = 5
@everywhere const mp = @SVector [2.4478382877847715e-06, 3.040432648022642e-06, 0.0009547919152112404, 0.0002858856727222417, 0.0]

@everywhere sun_mf_approx(r) = 0.5010980317741527*(2 + (-2 - 10.284885633196014*r*(2 + 10.284885633196014*r)) * exp(-10.284885633196014*r))

@everywhere function qdot!(dq,p,q,params,t)
	dq .= p
end

@everywhere function v3(q,i)
	return @inbounds SVector(q[1,i],q[2,i],q[3,i])
end

@everywhere function sv3(q,i,v)
	@inbounds begin
		q[1,i] += v[1]
		q[2,i] += v[2]
		q[3,i] += v[3]
	end
end

# Using MutableArray?
# Doing something more intelligent with test particles?
@everywhere function pdot!(dp,p,q,params,t)
	# Is this necessary?
	fill!(dp,0.)
	x3 = @SVector [norm(v3(q,i))^3 for i in 1:nplanets]
	for i in 1:nplanets
		xi = v3(q,i)
		#xi3 = norm(xi)^3
		@inbounds xi3 = x3[i]
		for j in i+1:nplanets
			xj = v3(q,j)
			# Recalc ...
			#xj3 = norm(xj)^3
			@inbounds xj3 = x3[j]
			dx = xi .- xj
			dx3 = norm(dx)^3
			@inbounds sv3(dp, i, -mp[j]*(dx/dx3 + xj/xj3))
			@inbounds sv3(dp, j, mp[i]*(dx/dx3 - xi/xi3))
			#@. dp[:,i] -= mp[j]*(dx/dx3 + xj/xj3)
			#@. dp[:,j] += mp[i]*(dx/dx3 - xi/xi3)
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
		#@. dp[:,i] -= pre*xi/xi3
	end
end


@everywhere function integ_diagnose(p,tmax,tgap,fn!)

	tspan = (0,1.5*tgap)

	t0 = p[1]
	tstep = p[2]
	x0p = zeros(3,5)
	v0p = zeros(3,5)

	for i in 1:5
		x0p[:,i] .= p[3 + 6*(i-1):3 + 6*(i-1) + 2]
		v0p[:,i] .= p[3 + 6*(i-1) + 3:3 + 6*(i-1) + 5]
	end

	while true
		prob = DynamicalODEProblem(pdot!,qdot!,v0p,x0p,tspan)
		integrator = init(prob,DPRKN12(),reltol=1e-14,abstol=1e-14,dt=tstep,maxiters=Inf)

		while true

			step!(integrator)
			err = check_error(integrator)
			if err != :Success
				println(err)
				return
			end

			fn!(integrator.t,integrator.u,integrator.tprev)

			if integrator.t > tgap
				t0::Float64 = t0 + integrator.t
				tstep::Float64 = integrator.t - integrator.tprev
				x0p .= integrator.u[2,:,:]
				v0p .= integrator.u[1,:,:]
				break
			end
		end

		if t0 > tmax 
			break
		end
	end
end

@everywhere function read_particle_csv_err_1(fname)
	open(fname) do f
		nlines = countlines(f)
		# Checking if we end in newline ...
		seekend(f)
		skip(f,-1)
		c = read(f,Char)
		if c != '\n'
			nlines -= 1
		end
		a = zeros(7,nlines)
		# back to beginning 
		seek(f,0)
		for i in 1:nlines
			line = readline(f)
			nums = split(line,",")
			t = parse(Float64,nums[1])
			xv = parse.(Float64,nums[27:32])
			a[1,i] = t
			a[2:7,i] .= xv
		end
		return a
	end
end

@everywhere function read_particle_csv_err_mem(fname)
	a = []
	lines = readlines(fname)
	for line in lines
		nums = split(line,",")
		if length(nums) != 32
			break
		end
		try
			t = parse(Float64,nums[1])
			xv = parse.(Float64,nums[27:32])
			push!(a,t)
			for f in xv
				push!(a,f)
			end
		catch err
			break
		end
	end
	reshape(a,7,:)
end

@everywhere function read_particle_csv_err(fname)
	a = []
	open(fname) do file
		while !eof(file)
			line = readline(file,keep=true)
			if last(line) != '\n'
				break
			end
			nums = split(line,",")
			if length(nums) != 32
				break
			end
			try
				t = parse(Float64,nums[1])
				xv = parse.(Float64,nums[27:32])
				push!(a,t)
				for f in xv
					push!(a,f)
				end
			catch err
				break
			end
		end
	end
	reshape(a,7,:)
end

@everywhere function read_particle_earth_csv_err_pos(f,p)
	seek(f,p)
	nlines = countlines(f)
	# Checking if we end in newline ...
	seekend(f)
	skip(f,-1)
	c = read(f,Char)
	if c != '\n'
		nlines -= 1
	end
	a = zeros(13,nlines)
	# back to beginning 
	seek(f,p)
	for i in 1:nlines
		line = readline(f)
		nums = split(line,",")
		t = parse(Float64,nums[1])
		xv = parse.(Float64,nums[27:32])
		xve = parse.(Float64,nums[9:14])
		a[1,i] = t
		a[2:7,i] .= xv
		a[8:13,i] .= xve
	end
	return a
end

@everywhere function read_particle_earth_csv_pos(f,p,pos_end)
	a = []
	seek(f,p)
	while position(f) != pos_end
		@assert !eof(f)
		line = readline(f)
		nums = split(line,",")
		t = parse(Float64,nums[1])
		xv = parse.(Float64,nums[27:32])
		xve = parse.(Float64,nums[9:14])
		push!(a,t)
		append!(a,xv)
		append!(a,xve)
	end
	reshape(a,13,:)
end

@everywhere function read_full_csv_pos(f,p,pos_end)
	a = []
	seek(f,p)
	while position(f) != pos_end
		@assert !eof(f)
		line = readline(f)
		nums = split(line,",")
		append!(a,parse.(Float64,nums))
	end
	reshape(a,32,:)
end

@everywhere function read_particle_csv_end(fname,pos_end)
	a = []
	open(fname) do file
		while position(file) != pos_end
			line = readline(file)
			nums = split(line,",")
			t = parse(Float64,nums[1])
			xv = parse.(Float64,nums[27:32])
			push!(a,t)
			append!(a,xv)
		end
	end
	reshape(a,7,:)
end

@everywhere function read_particle_earth_csv_err(fname)
	open(fname) do f
		nlines = countlines(f)
		# Checking if we end in newline ...
		seekend(f)
		skip(f,-1)
		c = read(f,Char)
		if c != '\n'
			nlines -= 1
		end
		a = zeros(13,nlines)
		# back to beginning 
		seek(f,0)
		for i in 1:nlines
			line = readline(f)
			nums = split(line,",")
			t = parse(Float64,nums[1])
			xv = parse.(Float64,nums[27:32])
			xve = parse.(Float64,nums[9:14])
			a[1,i] = t
			a[2:7,i] .= xv
			a[8:13,i] .= xve
		end
		return a
	end
end

function read_full_csv_err(fname)
	open(fname) do f
		nlines = countlines(f)
		# Checking if we end in newline ...
		seekend(f)
		skip(f,-1)
		c = read(f,Char)
		if c != '\n'
			nlines -= 1
		end
		a = zeros(32,nlines)
		# back to beginning 
		seek(f,0)
		for i in 1:nlines
			line = readline(f)
			nums = split(line,",")
			a[:,i] .= parse.(Float64,nums)
		end
		return a
	end
end

@everywhere function find_time(t,f,pos_end)
	# First line
	seek(f,0)
	s = readuntil(f,',')
	p0 = 0
	t0 = parse(Float64,s)
	if t < t0
		return nothing
	end
	# Then, check whether we end in newline
	# FIXME - shouldn't need to do this any more ...
	seek(f,pos_end)
	p1 = position(f)
	skipb = max(-4000,-(p1-1))
	skip(f,skipb)
	readline(f)
	# From here, going through lines ...
	p1,t1 = p0,t0
	while true
		p11 = position(f)
		#@assert !eof(f)
		@assert p11 <= pos_end
		if p11 == pos_end
			break
		end
		l = readline(f,keep=true)
		if last(l) != '\n'
			break
		end
		# Otherwise, accept as valid ...
		s = split(l,',')
		p1 = p11
		t1 = parse(Float64,s[1])
		if t1 > t
			break
		end
	end
	if t > t1
		return nothing
	end
	while true
		p2a = floor(Int,(p0+p1)/2)
		seek(f,p2a)
		skipchars(!=('\n'),f)
		skip(f,1)
		p2 = position(f)
		if p2 == p1
			return (p1,t1)
		elseif p2 == p0
			return (p0,t0)
		end
		s = readuntil(f,',')
		t2 = parse(Float64,s)
		if t < t2
			t1 = t2
			p1 = p2
		else
			t0 = t2
			p0 = p2
		end
	end
end

# So, finding particular time through binary search
@everywhere function find_time(t,f)
	# First line
	seek(f,0)
	s = readuntil(f,',')
	p0 = 0
	t0 = parse(Float64,s)
	if t < t0
		return nothing
	end
	# Then, check whether we end in newline
	seekend(f)
	p1 = position(f)
	skipb = max(-4000,-(p1-1))
	skip(f,skipb)
	readline(f)
	# From here, going through lines ...
	p1,t1 = p0,t0
	while !eof(f)
		p11 = position(f)
		l = readline(f,keep=true)
		if last(l) != '\n'
			break
		end
		# Otherwise, accept as valid ...
		s = split(l,',')
		p1 = p11
		t1 = parse(Float64,s[1])
		if t1 > t
			break
		end
	end
	if t > t1
		return nothing
	end
	while true
		p2a = floor(Int,(p0+p1)/2)
		seek(f,p2a)
		skipchars(!=('\n'),f)
		skip(f,1)
		p2 = position(f)
		if p2 == p1
			return (p1,t1)
		elseif p2 == p0
			return (p0,t0)
		end
		s = readuntil(f,',')
		t2 = parse(Float64,s)
		if t < t2
			t1 = t2
			p1 = p2
		else
			t0 = t2
			p0 = p2
		end
	end
end

# Firstly, looking for intersections with earth band

function xvp(a)
	(SVector(a[27],a[28],a[29]),SVector(a[30],a[31],a[32]))
end

@everywhere function xvi(a,i)
	b = 2 + 6*(i-1)
	(SVector(a[b+1],a[b+2],a[b+3]),SVector(a[b+4],a[b+5],a[b+6]))
end

# Assuming that it's just particle and Earth
@everywhere function shell_hits_rot(x,v,xe,ve)
	#x,v = xvp(ap)
	#xe, ve = xvi(ap,2)
	# Firstly, Earth orbits stuff ...
	# TODO - does radius change much?
	#earth_e = 0.5*dot(ve,ve) - 1/norm(xe)
	#earth_a = -0.5/earth_e
	earth_l = cross(xe,ve)
	earth_ev = cross(ve,earth_l) - xe/norm(xe)
	xp = earth_ev/norm(earth_ev)
	zp = earth_l/norm(earth_l)
	yp = cross(xp,zp)
	rmat = @SMatrix [xp[1] xp[2] xp[3];
					 yp[1] yp[2] yp[3];
					 zp[1] zp[2] zp[3]]
	r1xv(rmat * x, rmat * v)
end

@everywhere function shell_hits_rot(a)
	x = SVector(a[2],a[3],a[4]); v = SVector(a[5],a[6],a[7])
	xe = SVector(a[8],a[9],a[10]); ve = SVector(a[11],a[12],a[13])
	#x,v = xvp(ap)
	#xe, ve = xvi(ap,2)
	# Firstly, Earth orbits stuff ...
	# TODO - does radius change much?
	#earth_e = 0.5*dot(ve,ve) - 1/norm(xe)
	#earth_a = -0.5/earth_e
	earth_l = cross(xe,ve)
	earth_ev = cross(ve,earth_l) - xe/norm(xe)
	xp = earth_ev/norm(earth_ev)
	zp = earth_l/norm(earth_l)
	yp = cross(xp,zp)
	rmat = @SMatrix [xp[1] xp[2] xp[3];
					 yp[1] yp[2] yp[3];
					 zp[1] zp[2] zp[3]]
	r1xv(rmat * x, rmat * v)
end

@everywhere function hit_acc_file_tmin(fname,dcuts,tmin,pos_end)
	open(fname,"r") do f
		tf = find_time(tmin,f,pos_end)
		if tf == nothing
			return nothing
		end
		a = read_particle_earth_csv_pos(f,tf[1],pos_end)
		return hit_acc_arr(a,dcuts)
	end
end

@everywhere function hit_acc_file_tmin_integ(fname,dcuts,tmin,pos_end)
	open(fname,"r") do f
		tf = find_time(tmin,f,pos_end)
		if tf == nothing
			return nothing
		end
		a = read_full_csv_pos(f,tf[1],pos_end)
		#return hit_acc_integ(a,dcuts)
		return hit_acc_integ_all(a,dcuts)
		#return hit_acc_integ_all_shell(a,dcuts)
	end
end

@everywhere function hit_acc_file_tmin_integ_r(fname,bz,dc,tmin,pos_end)
	open(fname,"r") do f
		tf = find_time(tmin,f,pos_end)
		if tf == nothing
			return nothing
		end
		a = read_full_csv_pos(f,tf[1],pos_end)
		return hit_acc_integ_r(a,bz,dc)
		#return acc_integ_r(a,bz)
	end
end

@everywhere function hit_acc_file_tmin_integ_z0r(fname,bz,dc,tmin,pos_end)
	open(fname,"r") do f
		tf = find_time(tmin,f,pos_end)
		if tf == nothing
			return nothing
		end
		a = read_full_csv_pos(f,tf[1],pos_end)
		return acc_ecliptic_r(a,bz)
	end
end

@everywhere function hit_acc_file_tmin_integ_M(fname,bz,dc,tmin,pos_end)
	open(fname,"r") do f
		tf = find_time(tmin,f,pos_end)
		if tf == nothing
			return nothing
		end
		a = read_full_csv_pos(f,tf[1],pos_end)
		#return hit_acc_integ_M(a,bz,dc)
		return hit_acc_integ_M_all(a,bz,dc)
	end
end

@everywhere function hit_acc_file(f,dcuts)
	a = read_particle_earth_csv_err(f)
	hit_acc_arr(a,dcuts)
end

@everywhere function a_xv(x,v)
	e = 0.5*dot(v,v) - 1/norm(x)
	-0.5/e
end

# FIXME - handle initial locations in Sun
@everywhere function hit_acc_arr(a,dcuts)
	nt = size(a)[2]
	#dcuts = (0.1,0.01,0.005,0.002,0.001)
	acc = zeros(length(dcuts))
	t0  = a[1,1]
	for i in 2:nt
		ap = a[:,i]
		x = SVector(ap[2],ap[3],ap[4]); v = SVector(ap[5],ap[6],ap[7])
		if norm(x) < rsun
			continue
		end
		dt = ap[1] - t0
		r1l = shell_hits_rot(ap)
		if r1l == nothing
			continue
		end
		# So, we have hits ... 
		av = a_xv(x,v)
		torbit = 2pi * av^(3/2)
		for ir in 1:2
			x1,v1 = r1l[ir]
			xz = x1[3]
			vr = dot(v1,x1)/norm(x1)
			for (i,dc) in enumerate(dcuts)
				if abs(xz) < dc
					# Accumulate into ...
					acc[i] += dt/(abs(vr) * torbit)
				end
			end
		end
		t0 = ap[1]
	end
	acc
end

@everywhere function bidx(bz,x0,x1)
	u -> clamp(floor(Int,bz*(u-x0)/(x1-x0)),1,bz)
end

# Accumulate into velocity array ...
# TODO - doing more proper version using Earth position?
@everywhere function hit_vel_acc(a,dc)
	nt = size(a)[2]
	t0  = a[1,1]
	bz = 512
	ifn = bidx(bz,-1.5,1.5)
	ifn1 = bidx(bz,0.,1.5)
	acc = zeros(bz,bz)
	for i in 2:nt
		ap = a[:,i]
		x = SVector(ap[2],ap[3],ap[4]); v = SVector(ap[5],ap[6],ap[7])
		if norm(x) < rsun
			continue
		end
		dt = ap[1] - t0
		r1l = shell_hits_rot(ap)
		if r1l == nothing
			continue
		end
		# So, we have hits ... 
		av = a_xv(x,v)
		torbit = 2pi * av^(3/2)
		for ir in 1:2
			x1,v1 = r1l[ir]
			xz = x1[3]
			if abs(xz) < dc
				# Accumulate into ...
				vth = dot(v1,cross(SVector(0,0,1.),x1))
				vp = sqrt(dot(v1,v1) - vth^2)
				vr = dot(v1,x1)/norm(x1)
				# So, what's the correct quantity to accumulate here?
				acc[ifn(vth),ifn1(vp)] += dt/(abs(vr)*torbit)
				#acc[i] += dt/(abs(vr) * torbit)
			end
		end
		t0 = ap[1]
	end
	acc
end

#################################################

# So, let's do a thing where, if we get close enough to Earth, then we 
# do an actual simulation and see how close we get ... stepping through, analysing ...

@everywhere function hit_acc_integ(a,dcuts)
	nt = size(a)[2]
	t0  = a[1,1]
	acc = zeros(length(dcuts))
	for i in 2:nt
		ap = a[:,i]
		x, v = xvi(ap,5)
		if norm(x) < rsun
			continue
		end
		dt = ap[1] - t0
		t0 = ap[1]
		xe0, ve0 = xvi(ap,2)
		r1l = shell_hits_rot(x,v,xe0,ve0)
		if r1l == nothing
			continue
		end
		z1 = abs(r1l[1][1][3])
		z2 = abs(r1l[2][1][3])
		#if z1 > 0.02 && z2 > 0.02
		if z1 > 0.1 && z2 > 0.1
		#if z1 > 0.04 && z2 > 0.04
			continue
		end
		#acc[1] += 1.
		ap[1] = 0.
		function accfn(t,u,tp)
			xp = SVector(u[2,1,5],u[2,2,5],u[2,3,5])
			xe = SVector(u[2,1,2],u[2,2,2],u[2,3,2])
			d = norm(xp-xe)
			if d < dcuts[1]
				for (j,dc) in enumerate(dcuts)
					if d < dc
						acc[j] += (t-tp)
					end
				end
			end
		end
		integ_diagnose(ap,dt,1e3*year,accfn)
	end
	acc
end

@everywhere function hit_acc_integ_all(a,dcuts)
	nt = size(a)[2]
	t0  = a[1,1]
	acc = zeros(length(dcuts))
	for i in 2:nt
		ap = a[:,i]
		dt = ap[1] - t0
		t0 = ap[1]
		ap[1] = 0.
		function accfn(t,u,tp)
			xp = SVector(u[2,1,5],u[2,2,5],u[2,3,5])
			xe = SVector(u[2,1,2],u[2,2,2],u[2,3,2])
			d = norm(xp-xe)
			if d < dcuts[1]
				for (j,dc) in enumerate(dcuts)
					if d < dc
						acc[j] += (t-tp)
					end
				end
			end
		end
		integ_diagnose(ap,dt,1e3*year,accfn)
	end
	acc
end

@everywhere function hit_acc_integ_all_shell(a,dcuts)
	nt = size(a)[2]
	t0  = a[1,1]
	acc = zeros(length(dcuts))
	for i in 2:nt
		ap = a[:,i]
		dt = ap[1] - t0
		t0 = ap[1]
		ap[1] = 0.
		function accfn(t,u,tp)
			xp = SVector(u[2,1,5],u[2,2,5],u[2,3,5])
			d = abs(1 - norm(xp))
			if d < dcuts[1]
				for (j,dc) in enumerate(dcuts)
					if d < dc
						acc[j] += (t-tp)
					end
				end
			end
		end
		integ_diagnose(ap,dt,1e3*year,accfn)
	end
	acc
end

@everywhere function hit_acc_integ_r(a,bz,dc)
	nt = size(a)[2]
	t0  = a[1,1]
	epe = 0.016697099474039835
	ae = 1.0000263972072314
	ifn = bidx(bz,ae*(1-3*epe),ae*(1+3*epe))
	acc = zeros(bz)
	for i in 2:nt
		ap = a[:,i]
		x, v = xvi(ap,5)
		if norm(x) < rsun
			continue
		end
		dt = ap[1] - t0
		t0 = ap[1]
		xe0, ve0 = xvi(ap,2)
		r1l = shell_hits_rot(x,v,xe0,ve0)
		if r1l == nothing
			continue
		end
		z1 = abs(r1l[1][1][3])
		z2 = abs(r1l[2][1][3])
		if z1 > 0.05 && z2 > 0.05
			continue
		end
		function accfn(t,u,tp)
			xp = SVector(u[2,1,5],u[2,2,5],u[2,3,5])
			xe = SVector(u[2,1,2],u[2,2,2],u[2,3,2])
			re = norm(xe)
			d = norm(xp-xe)
			if d < dc
				acc[ifn(re)] += (t-tp)
			end
		end
		integ_diagnose(ap,dt,1e3*year,accfn)
	end
	acc
end

@everywhere rotz(t) = @SMatrix [cos(t) sin(t) 1.;
								-sin(t) cos(t) 0.;
								0. 0. 1.]

# Next, accumulating into bins depending on the distance of Earth from Sun ...
@everywhere function hit_acc_integ_M(a,bz,dc)
	nt = size(a)[2]
	t0  = a[1,1]
	ifn = bidx(bz,0,2pi)
	acc = zeros(bz)
	for i in 2:nt
		ap = a[:,i]
		x, v = xvi(ap,5)
		if norm(x) < rsun
			continue
		end
		dt = ap[1] - t0
		t0 = ap[1]
		xe0, ve0 = xvi(ap,2)
		r1l = shell_hits_rot(x,v,xe0,ve0)
		if r1l == nothing
			continue
		end
		z1 = abs(r1l[1][1][3])
		z2 = abs(r1l[2][1][3])
		if z1 > 0.05 && z2 > 0.05
			continue
		end
		function accfn(t,u,tp)
			xp = SVector(u[2,1,5],u[2,2,5],u[2,3,5])
			xe = SVector(u[2,1,2],u[2,2,2],u[2,3,2])
			ve = SVector(u[1,1,2],u[1,2,2],u[1,3,2])
			# rotation hack
			#rm = rotz(70 * pi/180)
			#xe = rm * xe
			#ve = rm * ve
			#
			d = norm(xp-xe)
			if d < dc
				M = mean_anomly(xe,ve)
				acc[ifn(M)] += (t-tp)
			end
		end
		integ_diagnose(ap,dt,1e3*year,accfn)
	end
	acc
end

@everywhere function hit_acc_integ_M_all(a,bz,dc)
	nt = size(a)[2]
	t0  = a[1,1]
	ifn = bidx(bz,0,2pi)
	acc = zeros(bz)
	for i in 2:nt
		ap = a[:,i]
		dt = ap[1] - t0
		t0 = ap[1]
		function accfn(t,u,tp)
			xp = SVector(u[2,1,5],u[2,2,5],u[2,3,5])
			xe = SVector(u[2,1,2],u[2,2,2],u[2,3,2])
			ve = SVector(u[1,1,2],u[1,2,2],u[1,3,2])
			d = norm(xp-xe)
			if d < dc
				M = mean_anomly(xe,ve)
				acc[ifn(M)] += (t-tp)
			end
		end
		integ_diagnose(ap,dt,1e3*year,accfn)
	end
	acc
end
#################################################

@everywhere function aupdown(a)
	nt = size(a)[2]
	bz = 512
	aup = zeros(Int,bz,bz)
	adown = zeros(Int,bz,bz)
	ifnx = bidx(bz,-1.,1.)
	ifny = bidx(bz,-7.,1.)
	aprev = 0.
	for i in 1:nt
		ap = a[:,i]
		x = SVector(ap[2],ap[3],ap[4]); v = SVector(ap[5],ap[6],ap[7])
		if norm(x) < rsun
			# Fixing this ...
			continue
		end
		av = a_xv(x,v)
		if av < 0.
			continue
		end
		if i == 1
			aprev = av
			continue
		end
		if av > aprev
			da = av - aprev
			aup[ifnx(log10(av)),ifny(log10(da))] += 1
		else
			da = aprev - av
			adown[ifnx(log10(av)),ifny(log10(da))] += 1
		end
		aprev = av
	end
	(aup, adown)
end

@everywhere function hit_sun_acc_arr(a)
	nt = size(a)[2]
	acc = 0.
	for i in 2:nt
		dt = a[1,i] - a[1,i-1]
		x = SVector(a[2,i],a[3,i],a[4,i]); v = SVector(a[5,i],a[6,i],a[7,i])
		# So, finding stuff ... 
		r = norm(x)
		e = 0.5*dot(v,v) + sun_phi(r)
		l = cross(x,v)
		l2 = dot(l,l)
		if r < rsun || (-1 + sqrt(1 + 2*e*l2))/(2*e) < rsun
			av = -0.5/e
			if av < 0.
				continue
			end
			nhits = dt/(pi*(av^(3/2)))
			vr1 = sqrt(2*(e + 1/rsun - 0.5*l2/rsun^2))
			acc += nhits/vr1
		end
	end
	acc
end


@everywhere function acc_integ_r(a,bz)
	nt = size(a)[2]
	acc = zeros(bz)
	epe = 0.016697099474039835
	ae = 1.0000263972072314
	# Check whether this is accurate ...
	rpts = range(ae*(1-3*epe),ae*(1+3*epe),length=bz)
	for i in 2:nt
		dt = a[1,i] - a[1,i-1]
		ap = a[:,i]
		xe, ve = xvi(ap,2)
		ee = (0.5*dot(ve,ve) - 1/norm(xe))
		ae = -0.5/ee
		le = norm(cross(xe,ve))
		ecc = sqrt(1 + 2*ee*le^2)
		torbit = 2pi*ae^(3/2)
		norbit = dt/torbit
		for j in 1:bz
			r = rpts[j]
			if r < ae*(1 + ecc) && r > ae*(1-ecc)
				vr = sqrt(2*(ee - 0.5*le^2/r^2 + 1/r))
				acc[j] += norbit / vr
			end
		end
	end
	acc
end


@everywhere function z0xv(x,v)
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
	# So, solving ...
	az = av[3]; bz = bv[3]
	tx =(az^2*ecc - sqrt(az^2*bz^2 + bz^4 - az^2*bz^2*ecc^2))/(az^2 + bz^2) 
	ty = (az*ecc - (az^3*ecc)/(az^2 + bz^2) + (
 az*sqrt(-bz^2*(-az^2 - bz^2 + az^2*ecc^2)))/(az^2 + bz^2))/bz
	t1 = atan(ty,tx)
	tx =(az^2*ecc + sqrt(az^2*bz^2 + bz^4 - az^2*bz^2*ecc^2))/(az^2 + bz^2) 
	ty = (-az*ecc + (az^3*ecc)/(az^2 + bz^2) + (
 az*sqrt(-bz^2*(-az^2 - bz^2 + az^2*ecc^2)))/(az^2 + bz^2))/bz
	t2 = atan(ty,tx)
	xv = zeros(3,2,2)
	for i in 1:2
		t = i == 1 ? t1 : t2
		cn = cos(t)
		sn = sin(t)
		dndt = 1/(a^1.5 * (1 - ecc*cn))
		xv[:,1,i] .= (-ecc + cn)*av + sn*bv
		xv[:,2,i] .= dndt*(-sn*av + cn*bv)
	end
	((xv[:,1,1],xv[:,2,1]),
	 (xv[:,1,2],xv[:,2,2]))
end

@everywhere function acc_ecliptic_r(a,bz)
	nt = size(a)[2]
	acc = zeros(bz)
	ifn = bidx(bz,0.4,4.)
	for i in 2:nt
		dt = a[1,i] - a[1,i-1]
		ap = a[:,i]
		x, v = xvi(ap,5)
		# Finding where we intersect the z=0 plane
		r1l = z0xv(x,v)
		if r1l == nothing
			continue
		end
		av = a_xv(x,v)
		torbit = 2pi * av^(3/2)
		for ir in 1:2
			x1,v1 = r1l[ir]
			r1 = norm(x1)
			vz = abs(v1[3]) 
			# Dealing with every small vz situation, then ...
			denom = vz*torbit
			if denom < 0.01
				denom = 0.01
			end
			# So, accumulate ...
			acc[ifn(r1)] += dt/denom
		end
	end
	acc
end

@everywhere function hit_sun_shell_acc_arr(a,rs)
	nt = size(a)[2]
	acc = 0.
	phis = sun_phi(rs)
	for i in 2:nt
		dt = a[1,i] - a[1,i-1]
		x = SVector(a[2,i],a[3,i],a[4,i]); v = SVector(a[5,i],a[6,i],a[7,i])
		# So, finding stuff ... 
		r = norm(x)
		e = 0.5*dot(v,v) + sun_phi(r)
		l = cross(x,v)
		l2 = dot(l,l)
		#
		es = 0.5*l2/rs^2 + phis
		if es < e
			av = -0.5/e
			if av < 0.
				continue
			end
			nhits = dt/(pi*(av^(3/2)))
			vr1 = sqrt(2*(e - es))
			acc += nhits/vr1
		end
	end
	acc
end

@everywhere function density_interp(l)
	x = sqrt(2) - l
	(0.032163859512083795*x - 0.2685163875228167*x^2 + 3.3339540527878104*x^3 - 
 13.593875661619942*x^4 + 38.92464292696844*x^5 - 46.7890088479447*x^6 + 
 25.82646211316403*x^7 - 5.523548048465026*x^8)
end

@everywhere function hit_sun_density_acc_arr(a)
	nt = size(a)[2]
	acc = 0.
	for i in 2:nt
		dt = a[1,i] - a[1,i-1]
		x = SVector(a[2,i],a[3,i],a[4,i]); v = SVector(a[5,i],a[6,i],a[7,i])
		# So, finding stuff ... 
		r = norm(x)
		e = 0.5*dot(v,v) + sun_phi(r)
		l = cross(x,v)
		l2 = dot(l,l)
		if r < rsun || (-1 + sqrt(1 + 2*e*l2))/(2*e) < rsun
			av = -0.5/e
			if av < 0.
				continue
			end
			nhits = dt/(pi*(av^(3/2)))
			# So, density stuff ...
			acc += nhits * density_interp(sqrt(l2/rsun))
		end
	end
	acc
end

# Generating plots here, then ... getting together ...

@everywhere function action_angle_vars(x,v)
	e = 0.5*dot(v,v) + sun_phi(norm(x))
	a = -0.5/e
	l = cross(x,v)
	if a < 0
		return (0.,0.,0.,0.,0.,0.)
	end
	J1 = sqrt(a)
	J2 = norm(l)
	J3 = l[3]
	# Angle variables, then ...
	ev = cross(v,l) - x/norm(x)
	zv = [0., 0., 1.]
	# Vector pointing towards ascending node ... check this ...
	n = cross(zv, l)
	omega = acos(dot(n,ev)/(norm(n)*norm(ev)))
	if ev[3] < 0
		omega = 2pi - omega
	end
	#
	Omega = acos(n[1]/norm(n))
	if n[2] < 0
		Omega = 2pi - Omega
	end
	# TODO - check this definition
	ecc = norm(ev)
	x1 = dot(x,ev)/ecc + ecc*a
	bv = cross(l,ev)
	y1 = dot(x,bv)/norm(bv)
	E = atan(y1/sqrt(1-ecc^2),x1)
	M = E - ecc * sin(E)
	(J1,M,J2,omega,J3,Omega)
end

@everywhere function mean_anomly(x,v)
	e = 0.5*dot(v,v) + sun_phi(norm(x))
	a = -0.5/e
	l = cross(x,v)
	if a < 0
		return (0.,0.,0.,0.,0.,0.)
	end
	ev = cross(v,l) - x/norm(x)
	ecc = norm(ev)
	x1 = dot(x,ev)/ecc + ecc*a
	bv = cross(l,ev)
	y1 = dot(x,bv)/norm(bv)
	E = atan(y1/sqrt(1-ecc^2),x1)
	M = E - ecc * sin(E)
	mod2pi(M)
end

# Do we want to do Jupiter thing? Coming back to ...

@everywhere function aa_bitmaps(a)
	nt = size(a)[2]
	bz = 512
    jm = 1.5
    fns = (bidx(bz,0.,jm),
           bidx(bz,0.,jm),
           bidx(bz,-jm,jm),
           bidx(bz,0.,2pi),
           bidx(bz,0.,2pi),
           bidx(bz,0.,2pi))
	pairs = ((1,2),(1,3),(1,4),(1,5),(1,6),
			 (2,3),(2,4),(2,5),(2,6),
			 (3,4),(3,5),(3,6),
			 (4,5),(4,6),
			 (5,6))
	np = length(pairs)
	aa = zeros(np,bz,bz)
	for i in 2:nt
		dt = a[1,i] - a[1,i-1]
		ap = a[:,i]
		x = SVector(ap[2],ap[3],ap[4]); v = SVector(ap[5],ap[6],ap[7])
		if norm(x) < rsun
			continue
		end
		aavars = action_angle_vars(x,v)
		for j in 1:np
			p1,p2 = pairs[j]
			xv = aavars[p1]
			yv = aavars[p2]
			aa[j,fns[p1](xv),fns[p2](yv)] += dt
		end
	end
	aa
end


# TODO - is this the best way to do this?
@everywhere function start_and_end_params_err(file)
	open(file) do f
		# Read first line
		l = readline(f)
		nums = split(l,",")
		as = parse.(Float64,nums)
		#
		seekend(f)
		pos1 = position(f)
		skip(f,-1)
		c = read(f,Char)
		#
		skipb = max(-2000,-(pos1-1))
		skip(f,skipb)
		l = readline(f)
		while !eof(f)
			l1 = readline(f)
			if eof(f) && c != '\n'
				break
			end
			l = l1
		end
		nums = split(l,",")
		ae = parse.(Float64,nums)
		return (as,ae)
	end
end

#=
function run_script_1()
	of0 = OutputFiles(1:256,"outsun")
	of = collect(of0)
	nf = length(of)
	#
	dcuts = (0.1,0.01,0.005,0.002,0.001)
	saccs = SharedArray(zeros(length(dcuts),nf))
	@sync @distributed for i in 1:nf
		f = of[i]
		saccs[:,i] .= hit_acc_file(f,dcuts)
	end
	accs = sdata(saccs)
	@save "accs1.jld2" dcuts accs
end

function run_script_2()
	of0 = OutputFiles(1:256,"outsun")
	of = collect(of0)
	nf = length(of)
	#
	dcuts = (0.1,0.01,0.005,0.002,0.001)
	tmin = (4.5e9 - 1e7)*year
	saccs = SharedArray(zeros(length(dcuts),nf))
	@sync @distributed for i in 1:nf
		f = of[i]
		hacc = hit_acc_file_tmin(f,dcuts,tmin)
		if hacc != nothing
			saccs[:,i] .= hacc
		end
	end
	accs = sdata(saccs)
	@save "accs3.jld2" dcuts accs
end
=#

function run_script_2b()
	@load "snapshot_sun_1.jld2" of pos_end tnow seps
	nf = length(of)
	dcuts = (0.1,0.01,0.005,0.002,0.001)
	tmin = (4.5e9 - 1e7)*year
	saccs = SharedArray(zeros(length(dcuts),nf))
	@sync @distributed for i in 1:nf
		f = of[i]
		tmi = max(tmin,seps[1,1,i])
		if tmi < seps[1,2,i]
			hacc = hit_acc_file_tmin(f,dcuts,tmi,pos_end[i])
			if hacc != nothing
				saccs[:,i] .= hacc
			end
		end
	end
	accs = sdata(saccs)
	@save "sun_acc_1e7.jld2" tnow accs
end

function run_script_2c()
	#@load "snapshot_sun_3.jld2" of pos_end tnow seps fnums
	@load "snapshot_sun_5.jld2" of pos_end tnow seps fnums
	nf = length(of)
	dcuts = (0.1,0.05,0.02,0.01,0.005,0.002,0.001)
	tmin = (4.5e9 - 1e6)*year
	saccs = SharedArray(zeros(length(dcuts),nf))
	@sync @distributed for i in 1:nf
		if fnums[i] <= 256
			f = of[i]
			tmi = max(tmin,seps[1,1,i])
			if tmi < seps[1,2,i]
				hacc = hit_acc_file_tmin_integ(f,dcuts,tmi,pos_end[i])
				#hacc = hit_acc_file_tmin(f,dcuts,tmi,pos_end[i])
				if hacc != nothing
					saccs[:,i] .= hacc
				end
			end
		end
	end
	accs = sdata(saccs)
	#@save "sun_acc_5_1e8_shell.jld2" tnow accs
	#@save "sun_acc_5_1e7_integ_mod.jld2" tnow accs
	@save "sun_acc_5_1e6_integ_mod2.jld2" tnow accs
	#@save "sun_acc_3_1e7_integ_2.jld2" tnow accs
	#@save "sun_acc_3_1e7_shell.jld2" tnow accs
end

function run_script_2cc()
	#@load "snapshot_sun_3.jld2" of pos_end tnow seps fnums
	@load "snapshot_sun_5.jld2" of pos_end tnow seps fnums
	nf = length(of)
	#nf = 30
	dcuts = (0.1,0.05,0.02,0.01,0.005,0.002,0.001)
	tmin = (4.5e9 - 1e7)*year
	saccs = SharedArray(zeros(length(dcuts),nf))
	@sync @distributed for i in 1:nf
		if fnums[i] <= 256
			f = of[i]
			tmi = max(tmin,seps[1,1,i])
			if tmi < seps[1,2,i]
				hacc = hit_acc_file_tmin_integ(f,dcuts,tmi,pos_end[i])
				#hacc = hit_acc_file_tmin(f,dcuts,tmi,pos_end[i])
				if hacc != nothing
					saccs[:,i] .= hacc
				end
			end
		end
	end
	accs = sdata(saccs)
	@save "sun_acc_5_1e7_integ_full.jld2" tnow accs
end

function run_script_2d()
	@load "snapshot_sun_3.jld2" of pos_end tnow seps fnums
	nf = length(of)
	tmin = (4.5e9 - 1e8)*year
	bz = 30
	dc = 0.02
	saccs = SharedArray(zeros(bz,nf))
	@sync @distributed for i in 1:nf
		if fnums[i] <= 256
			f = of[i]
			tmi = max(tmin,seps[1,1,i])
			if tmi < seps[1,2,i]
				hacc = hit_acc_file_tmin_integ_r(f,bz,dc,tmi,pos_end[i])
				if hacc != nothing
					saccs[:,i] .= hacc
				end
			end
		end
	end
	accs = sdata(saccs)
	epe = 0.016697099474039835
	ae = 1.0000263972072314
	rpts = range(ae*(1-3*epe),ae*(1+3*epe),length=bz)
	@save "sun_acc_3_1e8_integ_r.jld2" tnow accs rpts
	#@save "sun_acc_3_1e7_integ_r_nohit.jld2" tnow accs rpts
end

function run_script_2ee()
	@load "snapshot_sun_5.jld2" of pos_end tnow seps fnums
	nf = length(of)
	tmin = (4.5e9 - 1e6)*year
	bz = 20
	dc = 0.02
	saccs = SharedArray(zeros(bz,nf))
	@sync @distributed for i in 1:nf
		if fnums[i] <= 256
			f = of[i]
			tmi = max(tmin,seps[1,1,i])
			if tmi < seps[1,2,i]
				hacc = hit_acc_file_tmin_integ_M(f,bz,dc,tmi,pos_end[i])
				if hacc != nothing
					saccs[:,i] .= hacc
				end
			end
		end
	end
	accs = sdata(saccs)
	@save "sun_acc_5_1e6_integ_M_all.jld2" tnow accs
end

function run_script_2e()
	@load "snapshot_sun_3.jld2" of pos_end tnow seps fnums
	nf = length(of)
	tmin = (4.5e9 - 1e8)*year
	bz = 20
	dc = 0.02
	saccs = SharedArray(zeros(bz,nf))
	@sync @distributed for i in 1:nf
		if fnums[i] <= 256
			f = of[i]
			tmi = max(tmin,seps[1,1,i])
			if tmi < seps[1,2,i]
				hacc = hit_acc_file_tmin_integ_M(f,bz,dc,tmi,pos_end[i])
				if hacc != nothing
					saccs[:,i] .= hacc
				end
			end
		end
	end
	accs = sdata(saccs)
	@save "sun_acc_3_1e8_integ_M_rot.jld2" tnow accs
end

function run_script_2f()
	@load "snapshot_sun_3.jld2" of pos_end tnow seps fnums
	nf = length(of)
	#nf = 10
	tmin = (4.5e9 - 1e9)*year
	bz = 100
	dc = 0.02
	saccs = SharedArray(zeros(bz,nf))
	@sync @distributed for i in 1:nf
		if fnums[i] <= 256
			f = of[i]
			tmi = max(tmin,seps[1,1,i])
			if tmi < seps[1,2,i]
				hacc = hit_acc_file_tmin_integ_z0r(f,bz,dc,tmi,pos_end[i])
				if hacc != nothing
					saccs[:,i] .= hacc
				end
			end
		end
	end
	accs = sdata(saccs)
	@save "sun_acc_3_1e8_integ_z0r.jld2" tnow accs
end

function run_script_3()
	of0 = OutputFiles(1:256,"outearth")
	of = collect(of0)
	nf = length(of)
	#
	saccs = SharedArray(zeros(nf))
	seps = SharedArray(zeros(7,2,nf))
	@sync @distributed for i in 1:nf
		f = of[i]
		a = read_particle_csv_err(f)
		saccs[i] = hit_sun_acc_arr(a)
		na = size(a)[2]
		if na != 0
			seps[:,1,i] .= a[:,1]
			seps[:,2,i] .= a[:,na]
		end
	end
	accs = sdata(saccs)
	eps = sdata(seps)
	@save "earth_accs1.jld2" of accs eps
end

function run_script_3b()
	@load "snapshot_earth_1.jld2" of pos_end tnow 
	nf = length(of)
	#
	saccs = SharedArray(zeros(nf))
	@sync @distributed for i in 1:nf
		f = of[i]
		a = read_particle_csv_end(f,pos_end[i])
		saccs[i] = hit_sun_shell_acc_arr(a,rsun)
	end
	accs = sdata(saccs)
	@save "earth_accs_1.jld2" tnow accs 
end

# For the a thing ... doing this for backwards runs, I guess?
# Probably a more illustrative sample space ... fine, doing that ...
# Accumulating into a things ...
function run_script_5()
	@load "snapshot_earth_2.jld2" of pos_end tnow 
	nf = length(of)
	rfn = ((a1,b1),(a2,b2)) -> (a1 .+ a2, b1 .+ b2)
	(aup,adown) = @distributed rfn for i in 1:nf
		f = of[i]
		a = read_particle_csv_end(f,pos_end[i])
		aupdown(a)
	end
	@save "earth_aupdown_2.jld2" aup adown
end

function run_script_6()
	@load "snapshot_sun_2.jld2" of pos_end tnow 
	nf = length(of)
	rfn = (b1,b2) -> b1 .+ b2
	aa = @distributed rfn for i in 1:nf
		f = of[i]
		a = read_particle_csv_end(f,pos_end[i])
		aa_bitmaps(a)
	end
	@save "sun_aagrid_2.jld2" aa
end

function run_script_4()
	@load "snapshot_sun_1.jld2" of pos_end tnow seps
	nf = length(of)
	dcut = 0.01
	tmin = (4.5e9 - 1e7)*year
	rfn = (b1,b2) -> b1 .+ b2
	bz = 512
	b = @distributed rfn for i in 1:nf
		f = of[i]
		tmi = max(tmin,seps[1,1,i])
		open(fname,"r") do f
			tf = find_time(tmi,f,pos_end[i])
			if tf == nothing
				return zeros(bz,bz)
			end
			a = read_particle_earth_csv_pos(f,tf[1],pos_end[i])
			hit_vel_acc(a,dcut)
		end
	end
	#@save "vel_acc_1e7_0p01.jld2" tnow b
	#=
	bmax = maximum(b)
	bout = ones(3,bz,bz)
	for i in CartesianIndices(b)
		bout[:,i] .-= b[i]/bmax
	end
	save("fig110.png",colorview(RGB,bout))
	=#
	bmax = maximum(b)
	# So, doing log thing ...
	bmin = bmax
	for i in eachindex(b)
		if b[i] != 0. && b[i] < bmin
			bmin = b[i]
		end
	end
	#
	b2 = (log.(b) .- log(bmin))/(log(bmax)-log(bmin))
	replace!(b2,-Inf=>0.)
	#
	bout = ones(3,bz,bz)
	# Now, doing the basic log thing, then ...
	for i in CartesianIndices(b)
		#bout[:,i] .-= b[i]/bmax
		bout[:,i] .-= b2[i]
	end
	save("fig114.png",colorview(RGB,bout))
end

#run_script_2c()
run_script_2cc()
#run_script_2ee()
#run_script_2d()
#run_script_2e()
#run_script_2f()
