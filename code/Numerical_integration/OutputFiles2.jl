# formats output filenames


struct OutputFiles
	irange
	dirbase
	dir_its
end

function OutputFiles(irange,dirbase)
	nd = length(irange)
	of = OutputFiles(irange,dirbase,zeros(Int,nd))
	# Filling in ...
	for (j,i) in enumerate(irange)
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
		of.dir_its[j] = ilast+1
	end
	of
end

function Base.iterate(of::OutputFiles)
	nd = of.irange[1]
	outstr = "$(of.dirbase)_$nd/output_$nd.csv"
	i = 1
	di = 1
	if di >= of.dir_its[i]
		i += 1
		di = 1
	else
		di += 1
	end
	return (outstr,(i,di))
end

function Base.iterate(of::OutputFiles,state)
	i,di = state
	if i > length(of.irange)
		return nothing
	end
	nd = of.irange[i]
	outstr = "$(of.dirbase)_$nd/output_$(1000*(di-1) + nd).csv"
	if di >= of.dir_its[i]
		i += 1
		di = 1
	else
		di += 1
	end
	return (outstr,(i,di))
end

Base.length(of::OutputFiles) = sum(of.dir_its)

function of_first(of::OutputFiles)
	OutputFiles(of.irange,of.dirbase,ones(Int,length(of.irange)))
end

#function of_finished(of::OutputFiles)
	#OutputFiles(of.ndirs,of.dirbase,max.(0,of.dir_its .- 1))
#end

# Zeroing out some of them ...

