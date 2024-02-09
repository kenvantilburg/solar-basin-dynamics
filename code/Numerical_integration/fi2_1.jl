include("SolarBasin2.jl")

#ENV["GKSwstype"]="nul"

using Plots
using Images
using Distributed

#upscale = 8 #8x upscaling in resolution
#fntsm = Plots.font("sans-serif", pointsize=round(10.0*upscale))
#fntlg = Plots.font("sans-serif", pointsize=round(14.0*upscale))
#default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
#default(size=(800*upscale,600*upscale)) #Plot canvas size

# Iterating through files ...

struct OutputFiles
	ndirs
	dirbase
	dir_its
end

function OutputFiles(ndirs,dirbase)
	of = OutputFiles(ndirs,dirbase,zeros(Int,ndirs))
	# Filling in ...
	for i in 1:ndirs
		fname = "$(dirbase)_$i/its"
		if !isfile(fname)
			break
		end
		ilast = 0
		open(fname,"r") do f
			for line in readlines(f)
				ilast = parse(Int,line)
			end
		end
		of.dir_its[i] = ilast+1
	end
	of
end

function Base.iterate(of::OutputFiles,state=(1,1))
	nd = state[1]
	di = state[2]
	if nd > of.ndirs
		return nothing
	end
	outstr = "$(of.dirbase)_$nd/output_$(1000*(di-1) + nd).csv"
	if di == of.dir_its[nd]
		nd += 1
		di = 1
	else
		di += 1
	end
	return (outstr,(nd,di))
end

Base.length(of::OutputFiles) = sum(of.dir_its)

function of_first(of::OutputFiles)
	OutputFiles(of.ndirs,of.dirbase,ones(Int,of.ndirs))
end

function of_finished(of::OutputFiles)
	a = []
	for nd in 1:of.ndirs
		ni = of.dir_its[nd]
		for di in 1:ni-1
			push!(a,"$(of.dirbase)_$nd/output_$(1000*(di-1) + nd).csv")
		end
	end
	a
end

function how_many(outbase, irange)
	sum = 0
	for i in irange
		fname = "$(outbase)_$i/its"
		ilast = 0
		open(fname,"r") do f
			for line in readlines(f)
				ilast = parse(Int,line)
			end
		end
		sum += ilast+1
	end
	sum
end

@everywhere function read_integ_csv2_err_stride(fname,s)
	# Make array, then reshape? Fixing ... 
	a = []
	open(fname) do file
		i = -1
		while !eof(file)
			line = readline(file)
			i += 1
			if (i % s != 0) && !eof(file)
				continue
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

@everywhere function ap(a)
	x = a[2:4]
	v = a[5:7]
	e = e = 0.5*dot(v,v) + sun_phi(norm(x))
	-0.5/e
end

function alz(a)
	x = a[2:4]
	v = a[5:7]
	cross(x,v)[3]
end

function xva(a)
	#(@SVector a[2:4], @SVector a[5:7])
	(SVector(a[2],a[3],a[4]),SVector(a[5],a[6],a[7]))
end

function bitmap_acc_params!(b,p,fx,fy)
	n = size(p)[2]
	for i in 1:n-1
		dt = p[1,i+1] - p[1,i]
		dt = max(dt,0.)
		u = p[:,i]
		b[fx(u),fy(u)] += dt
	end
end

@everywhere function bidx(bz,x0,x1)
	u -> clamp(floor(Int,bz*(u-x0)/(x1-x0)),1,bz)
end

function plota_bitmap(of,figfile,s)
	bz = 2048
	b = zeros(bz,bz)
	for f in of
		a = read_integ_csv2_err_stride(f,s)
		bitmap_acc_params!(b,a,
						  bidx(bz,0.,5.) ∘ ap,
						  bidx(bz,0.,2e9*year) ∘ first)
	end
	xvals = range(0,2e9,length=bz)
	yvals = range(0,5.,length=bz)
	plot()
	heatmap(xvals,yvals,log.(b),
		legend=:none,
		xlabel = "t/yr", ylabel = "a/AU", size=(2048,1536))
	png(figfile)
end

cpos(x) = max(x,0.)

function alz_bitmap(of,figfile,s)
	bz = 2048
	b = zeros(bz,bz)
	for f in of
		a = read_integ_csv2_err_stride(f,s)
		bitmap_acc_params!(b,a,
						  bidx(bz,-2.,2.) ∘ alz,
						  bidx(bz,0.,3.) ∘ sqrt ∘ cpos ∘ ap)
	end
	xvals = range(0,3.,length=bz)
	yvals = range(-2.,2.,length=bz)
	plot()
	heatmap(xvals,yvals,log.(b),
		legend=:none,
		xlabel = "a^1/2", ylabel = "Lz", size=(2048,1536))
	png(figfile)
end

#=

How to do processes stuff?

@spawnat (@spawn should work?)
@everywhere
fetch

Ah, okay, Base.Threads.@spawn
vs
Distributed.@spawnat
...

We need to use Distributed ...

@distributed ...
Distributed memory for loop ...

How're we going to do the reduction thing?
Either shared memory, or reduction thing ... how does reduction work?

x = f(x,v[i])

@distributed for, "final reduction done on calling process" ... okay ...
We can probably do a recursive reduction thing, I guess ...
Maybe easier to manage more explicitly?

"Transform collection c by applying f to each element" ...

SharedArray stuff ... ?

Ah, right, worker does reduction on stuff *that it's been assigned to*
Okay, good, that all makes sense, implementing like that, I guess?

Setting number of process with
julia -p n

Setting large-ish number here, then ..

=#

function plota_bparr(of,figfile,s)
	bz = 1024
	l = length(of)
	ba = u -> bz + 1 - bidx(bz,0.,5.)(u)
	bt =  bidx(bz,0.,3e9*year)
	ac = u -> u < 0. ? 5. : u
	#
	rfn = (b1,b2) -> min.(b1,b2)
	b = @distributed rfn for i in 1:l
		f = of[i]
		bf = ones(3,bz,bz)
		a = read_integ_csv2_err_stride(f,s)
		n = size(a)[2]
		# Random color
		ci = @SVector rand(3)
		ci *= 0.5
		for i in 1:n-1
			a0 = ac(ap(a[:,i]))
			a1 = ac(ap(a[:,i+1]))
			p0 = (ba(a0),bt(a[1,i]))
			p1 = (ba(a1),bt(a[1,i+1]))
			bline!(bf,ci,p0,p1)
		end
		bf
	end
	#
	cr = SVector(0.8,0.8,0.8)
	xrule!(b,cr,bt(1e9*year))
	xrule!(b,cr,bt(2e9*year))
	yrule!(b,cr,ba(1.))
	yrule!(b,cr,ba(2.))
	yrule!(b,cr,ba(3.))
	yrule!(b,cr,ba(4.))
	#
	save(figfile,colorview(RGB,b))
end

function alz_bitmap_parr(of,figfile,s)
	bz = 2048
	b = zeros(bz,bz)
	l = length(of)
	bi = zeros(bz,bz,l)
	Threads.@threads for i in 1:l
		f = of[i]
		bf = zeros(bz,bz)
		a = read_integ_csv2_err_stride(f,s)
		bitmap_acc_params!(bf,a,
						  bidx(bz,-2.,2.) ∘ alz,
						  bidx(bz,0.,3.) ∘ sqrt ∘ cpos ∘ ap)
		bi[:,:,i] .= bf
	end
	# Reduce
	for i in 1:l
		b .+= @view bi[:,:,i]
	end
	xvals = range(0,3.,length=bz)
	yvals = range(-2.,2.,length=bz)
	plot()
	heatmap(xvals,yvals,log.(b),
		legend=:none,
		xlabel = "a^1/2", ylabel = "Lz", size=(2048,1536))
	png(figfile)
end

function earth_vel_bitmap_parr(of,figfile,s)
	bz = 2048
	b = zeros(bz,bz)
	l = length(of)
	bi = zeros(bz,bz,l)
	ifn = bidx(bz,-1.5,1.5)
	Threads.@threads for i in 1:l
		f = of[i]
		bf = zeros(bz,bz)
		a = read_integ_csv2_err_stride(f,s)
		n = size(a)[2]
		for i in 1:n-1
			x,v = xva(a[:,i])
			r1l = r1locations(x,v)
			if r1l == nothing
				continue
			end
			dt = a[1,i+1] - a[1,i]
			dt = max(dt,0.)
			# Velocity relative to thing ...
			for ir in 1:2
				x1,v1 = r1l[ir]
				if abs(x1[3]) > 0.05
					continue
				end
				# If hits close enough to ecliptic, accumulate into
				# bitmap ...
				vr = dot(v1,x1)
				vth = dot(v1,cross(SVector(0,0,1.),x1))
				bf[ifn(vth),ifn(vr)] += dt/abs(vr)
			end
		end
		bi[:,:,i] .= bf
	end
	# Reduce
	for i in 1:l
		b .+= @view bi[:,:,i]
	end
	xvals = range(-1.5,1.5,length=bz)
	plot()
	heatmap(xvals,xvals,log.(b),
		legend=:none,
		xlabel = "vr", ylabel = "vt", size=(2048,2048))
	png(figfile)
end

# Bresenham, no anti-aliasing
@everywhere function bline!(b,c,(x0,y0),(x1,y1))
	dx = abs(x1-x0)
	sx = x0 < x1 ? 1 : -1
	dy = -abs(y1-y0)
	sy = y0 < y1 ? 1 : -1
	err = dx + dy
	while true
		b[:,x0,y0] .= c
		if x0 == x1 && y0 == y1
			break
		end
		e2 = 2*err
		if e2 >= dy
			err +=  dy
			x0 += sx
		end
		if e2 <= dx
			err += dx
			y0 += sy
		end
	end
end

function yrule!(b,c,x)
	s = size(b)
	for i in 1:s[3]
		b[:,x,i] .= c
	end
end

function xrule!(b,c,y)
	s = size(b)
	for i in 1:s[2]
		b[:,i,y] .= c
	end
end

#=
function plota_bitmap_test(of,figfile,s)
	bz = 1024
	b = zeros(bz,bz)
	for f in of
		a = read_integ_csv2_err_stride(f,s)
		bitmap_acc_params!(b,a,
						  bidx(bz,0.,5.) ∘ ap,
						  bidx(bz,0.,2e9*year) ∘ first)
	end
	xvals = range(0,2e9,length=bz)
	yvals = range(0,5.,length=bz)
	heatmap(xvals,yvals,log.(b),legend=:none,xlabel = "t/yr", ylabel = "a/AU")
	#heatmap(log.(b),legend=:none,aspect_ratio=1.,xlabel = "t/yr", ylabel = "a/AU")
	png(figfile)
end
=#

function plota(of,figfile,s)
	plot()
	for f in of
		a = read_integ_csv2_err_stride(f,s)
		n = size(a)[2]
		plot!(a[1,:]/year,[clamp(ap(a[:,i]),0.,5.) for i in 1:n],label="",
			  xlabel="t/yr",ylabel="a/AU")
	end
	png(figfile)
end

function start_time(file)
	s = ""
	open(file,"r") do f
		l = readline(f)
		nums = split(l,",")
		s = nums[1]
	end
	return parse(Float64,s)
end

function start_times(of)
	[start_time(f) for f in of]
end

function start_and_end_params(file)
	open(file) do f
		# Read first line
		l = readline(f)
		nums = split(l,",")
		as = parse.(Float64,nums)
		pos0 = position(f)
		# Read last line
		seekend(f)
		pos1 = position(f)
		if pos0 == pos1
			return (as,nothing)
		end
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
		ae = parse.(Float64,nums)
		return (as,ae)
	end
end

# Firstly, plot of ejected versus other stuff ...
function neject_plot(of,figfile)
	# Firstly, list of start and end times ...
	np = length(of)
	ti = zeros(np)
	ej = zeros(Bool,np)
	fin = zeros(Bool,np)
	for (i,f) in enumerate(of)
		(p0,p1) = start_and_end_params(f)
		if p1 == nothing
			continue
		end
		ti[i] = p1[1] - p0[1]
		ej[i] = norm(p1[27:32]) > 10.
		fin[i] = p1[1] >= 4.5e9*year
	end
	# So, have things ..
	a1 = np*ones(Int,np)
	a2 = zeros(Int,np)
	a3 = zeros(Int,np)
	a4 = zeros(Int,np)

	perm = sortperm(ti)
	nej = np
	nfin = np
	for (i,p) in enumerate(perm)
		nej -= ej[p]
		nfin -= fin[p] + ej[p]
		a2[i] = nej
		a3[i] = nfin 
	end

	a4 = [i-(np-a3[i]) for i in 1:np]

	ya = ti[perm] / year
	
	plot(size=(2048,1536))
	plot!(ya,a1,legend=:none,xaxis="t/year",yaxis="nparticles",
		guidefontsize=18,xtickfontsize=18,ytickfontsize=18)
	plot!(ya,a2,legend=:none)
	plot!(ya,a3,legend=:none)
	plot!(ya,a4,legend=:none)
	
	png(figfile)
end

function run_script()
    of = OutputFiles(1,"outearth")
    #of = OutputFiles(256,"outearth")
    of1 = collect(of)
    plota_bparr(of1,"fig64.png",10)
end

run_script()
