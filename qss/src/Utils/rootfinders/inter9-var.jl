#using BenchmarkTools

@inline function smallest(args::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds m = args[1]
	for i in 2:N
		@inbounds arg = args[i]
		m = arg < m ? arg : m
	end
	m
end

@inline function horner(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
	end
	v
end

@inline function horner2(x::F, y::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	w = v
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

@inline function hornerd(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	d = zero(F)
	for i in 2:N
		@inbounds coeff = coeffs[i]
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

@inline function intervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	cohigh = colow
	if low < zero(F)
		for i in 2:N
			@inbounds coeff = coeffs[i]
			colow, cohigh = if colow > zero(F)
				muladd(low, cohigh, coeff), muladd(high, colow, coeff)
			elseif cohigh < zero(F)
				muladd(high, cohigh, coeff), muladd(low, colow, coeff)
			else
				muladd(low, cohigh, coeff), muladd(low, colow, coeff)
			end
		end
	else
		for i in 2:N
			@inbounds coeff = coeffs[i]
			colow, cohigh = if colow > zero(F)
				muladd(low, colow, coeff), muladd(high, cohigh, coeff)
			elseif cohigh < zero(F)
				muladd(high, colow, coeff), muladd(low, cohigh, coeff)
			else
				muladd(high, colow, coeff), muladd(high, cohigh, coeff)
			end
		end	
	end
    colow, cohigh
end

@inline function posintervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	cohigh = colow
	for i in 2:N
		@inbounds coeff = coeffs[i]
		colow, cohigh = if colow > zero(F)
			muladd(low, colow, coeff), muladd(high, cohigh, coeff)
		elseif cohigh < zero(F)
			muladd(high, colow, coeff), muladd(low, cohigh, coeff)
		else
			muladd(high, colow, coeff), muladd(high, cohigh, coeff)
		end
	end
	colow, cohigh
end

function smallestpositiverootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, doms::Ptr{NTuple{2,F}}) where {N, F <: AbstractFloat}
	if N == 1
		return typemax(F)
	elseif N == 2
		@inbounds ret = -coeffs[2] / coeffs[1]
		if ret < zero(F)
			return typemax(F)
		else
			return ret
		end
	end
    @inbounds _coeff1 = inv(coeffs[1])
    poly = ntuple(N) do i
        @inbounds _coeff1 * coeffs[i]
    end
    poly′ = ntuple(N - 1) do i
        @inbounds (N-i) * poly[i]
    end
	MM = smallest(poly)
	if MM > zero(F) return typemax(F) end
	domlow, domhigh = zero(F), one(F) - MM
	index = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom′low, codom′high = posintervalhorner(domlow, domhigh, poly′)
		#@show comid codom′low codom′high
		if codom′low < zero(F) < codom′high
			leftlow, lefthigh, rightlow, righthigh = if comid < zero(F)
				domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
			else
				domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
			end
			#@show leftlow lefthigh rightlow righthigh
			if leftlow < lefthigh
				if rightlow < righthigh
					index += 1
					unsafe_store!(doms, (rightlow, righthigh), index)
				end
				domlow, domhigh = leftlow, lefthigh
				continue
			elseif rightlow < righthigh
				domlow, domhigh = rightlow, righthigh
				continue
			end
		else
			codomlow, codomhigh = horner2(domlow, domhigh, poly)
			#@show domlow domhigh codomlow codomhigh
			if codomlow * codomhigh < zero(F)
				while true
					comid, comid′ = hornerd(mid, poly)
					delta = comid / comid′
					newmid = mid - delta
					if abs(delta) < 1.0e-8mid 
						return newmid
					elseif domlow < newmid < domhigh
						mid = newmid
					else
						if comid * codomlow < zero(F)
							domhigh, codomhigh = mid, comid
						else
							domlow, codomlow = mid, comid
						end
						mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
					end
				end
			end
		end
		if index == 0 break end
		domlow, domhigh = unsafe_load(doms, index)
		index -= 1
	end
	return typemax(F)
end

function allrealrootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, res::Ptr{F}, doms::Ptr{NTuple{2,F}}) where {N, F <: AbstractFloat}
	if N == 1
		return 0
	elseif N == 2
		println("fix setIndex Ptr error")
		@inbounds res[1] = -coeffs[2] / coeffs[1]
		return 1
	end
    @inbounds _coeff1 = inv(coeffs[1])
    poly = ntuple(N) do i
        @inbounds _coeff1 * coeffs[i]
    end
    poly′ = ntuple(N - 1) do i
        @inbounds (N-i) * poly[i]
    end
	mm = zero(F)
	MM = zero(F)
	s = -one(F)
	for i in 2:N
		@inbounds coeff = poly[i]
		if coeff < MM; MM = coeff end
		_coeff = s * coeff
		if _coeff < mm; mm = _coeff end
		s = -s
	end
	if mm == zero(F) && MM == zero(F) return 0 end
	index = 0
	domlow, domhigh = if mm < zero(F)
		if MM < zero(F)
			index = 1
			unsafe_store!(doms, (zero(F), one(F) - MM), 1)
		end
		mm - one(F), zero(F)
	else
		zero(F), one(F) - MM
	end
    counter = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom′low, codom′high = intervalhorner(domlow, domhigh, poly′)
		#@show mid comid codom′low codom′high
		if codom′low < zero(F) < codom′high
			leftlow, lefthigh, rightlow, righthigh = if comid < zero(F)
				domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
			else
				domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
			end
			#@show leftlow lefthigh rightlow righthigh
			if leftlow < lefthigh
				if rightlow < righthigh
					index += 1
					unsafe_store!(doms, (rightlow, righthigh), index)
				end
				domlow, domhigh = leftlow, lefthigh
				continue
			elseif rightlow < righthigh
				domlow, domhigh = rightlow, righthigh
				continue
			end
		else
			codomlow, codomhigh = horner2(domlow, domhigh, poly)
            #@show domlow domhigh codomlow codomhigh
			if codomlow * codomhigh < zero(F)
				while true
					comid, comid′ = hornerd(mid, poly)
					delta = comid / comid′
                    #@show mid comid comid′
					newmid = mid - delta
					if abs(delta) < 1.0e-8abs(mid) 
						counter += 1
                        unsafe_store!(res, newmid, counter)
                        break
					elseif domlow < newmid < domhigh
						mid = newmid
					else
						if comid * codomlow < zero(F)
							domhigh, codomhigh = mid, comid
						else
							domlow, codomlow = mid, comid
						end
						mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
					end
				end
			end
		end
		if index == 0 || counter == N-1 break end
		domlow, domhigh = unsafe_load(doms, index)
		index -= 1
	end
	return counter
end








function classicRoot(coeff::NTuple{3,Float64}) # credit goes to github.com/CIFASIS/qss-solver
    mpr=(-1.0,-1.0) #
	a=coeff[1];b=coeff[2];c=coeff[3]
    if a == 0.0 #|| (10000 * abs(a)) < abs(b)# coef3 is the coef of t^2
        if b != 0.0
			if -c / b>0.0
			 mpr = (-c / b,-1.0)
			end
		  end
       
    else 
       #double disc;
        disc = b * b - 4.0 * a * c#b^2-4ac
        if disc > 0.0 # no real roots
         
        
          #double sd, r1;
          sd = sqrt(disc);
		  r1 = (-b + sd) / (2.0 * a);
       
          r2 = (-b - sd) / (2.0 * a);

			mpr = (r1,r2)
		
		elseif disc == 0.0
			r1 = (-b ) / (2.0 * a);
			mpr = (r1,-1.0)
		end
        
    end
    return mpr
end


function quadRootv2(coeff::NTuple{3,Float64}) # 
	mpr=(-1.0,-1.0) #size 2 to use mpr[2] in quantizer
	a=coeff[1];b=coeff[2];c=coeff[3]
	if a == 0.0 #|| (10000 * abs(a)) < abs(b)# coef3 is the coef of t^2
		if b != 0.0
		  if -c / b>0.0
		   mpr = (-c / b,-1.0)
		  end
		end
	elseif b==0.0
		if -c/a>0
		mpr = (sqrt(-c / a),-1.0)
		end
	elseif c==0.0
		mpr=(-1.0,-b/a)
	else 
	   #double disc;
	   Δ = 1.0 - 4.0*c*a / (b*b)
		if Δ >0.0
			#= q = -0.5*(1.0+sign(b)*sqrt(Δ))*b
			r1 = q / a
		   
			r2=c / q =#
			sq=sqrt(Δ)
			r1=-0.5*(1.0+sq)*b/a
			r2=-0.5*(1.0-sq)*b/a
		 
			mpr = (r1,r2)
		elseif Δ ==0.0
			r1=-0.5*b/a
			mpr = (r1,r1-1e-12)
		end
	end
	return mpr
  end



  function mprv2(coeff::NTuple{3,Float64}) # 
	mpr=Inf
	a=coeff[1];b=coeff[2];c=coeff[3]
	if a == 0.0 #|| (10000 * abs(a)) < abs(b)# coef3 is the coef of t^2
		if b != 0.0
		  if -c / b>0.0
		   mpr = -c / b
		  end
		end
	elseif b==0.0
		if -c/a>0
		mpr = sqrt(-c / a)
		end
	elseif c==0.0
		if -b/a>0.0
			mpr = -b/a
		end
		
	else 
	   #double disc;
	   Δ = 1.0 - 4.0*c*a / (b*b)
		if Δ >0.0
			#= q = -0.5*(1.0+sign(b)*sqrt(Δ))*b
			r1 = q / a
		   
			r2=c / q =#
			sq=sqrt(Δ)
			r1=-0.5*(1.0+sq)*b/a
			if r1 > 0 
				mpr = r1;
			end
			r1=-0.5*(1.0-sq)*b/a

			if ((r1 > 0) && (r1 < mpr)) 
				mpr = r1;
			  end
		elseif Δ ==0.0
			r1=-0.5*b/a
			if r1 > 0 
				mpr = r1;
			end
		end
	end
	

	return mpr
  end


  
				#= coeffs2=NTuple{3,Float64}((14.691504647354595,-747452.6968034876,1.0e-6))
						function iter(res1::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},coeffs2::NTuple{3,Float64})
											#coeffs2=NTuple{3,Float64}((1.0, -2.0, -1.06))
										
										pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
										res1 = pointer(Vector{Float64}(undef, 2))
										unsafe_store!(res1, -1.0, 1);unsafe_store!(res1, -1.0, 2)
										allrealrootintervalnewtonregulafalsi(coeffs2,res1,pp)
										resTup=(unsafe_load(res1,1),unsafe_load(res1,2))
										#@show resTup
										resfilterd=filter((x) -> x >0.0 , resTup)
										
									#	display(resfilterd)
						end
						function  anal1(coeffs2::NTuple{3,Float64})
							#coeffs2=NTuple{3,Float64}((1.0, -2.0, -1.06))
						    classicRoot(coeffs2)
						end
						function  anal2(coeffs2::NTuple{3,Float64})
							#coeffs2=NTuple{3,Float64}((1.0, -2.0, -1.06))
						    quadRootv2(coeffs2)
						end

						pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
						res1 = pointer(Vector{Float64}(undef, 2))
						@show iter(res1,pp,coeffs2)
						@show anal1(coeffs2)
						@show anal2(coeffs2)  =#

#= @btime iter(res1,pp)
@btime anal1()
@btime anal2() =#


#= coeffs=NTuple{3,Float64}((29160.956861496, 67.56376290717117, 0.014452560693316113))
#coeffs2=NTuple{3,Float64}((201011.0777843986, 106.4755863000737, 0.014452560693316113))
coeffs2=NTuple{4,Float64}((1.021243315770821e8, -3914.116824448214, -0.040323840000000166,1e-6))
pp=pointer(Vector{NTuple{2,Float64}}(undef, 11))
res1 = pointer(Vector{Float64}(undef, 3))
res2 = pointer(Vector{Float64}(undef, 3)) =#
#= count=allrealrootintervalnewtonregulafalsi(coeffs2,res2,pp)
@show count
resTup2=(unsafe_load(res2,1),unsafe_load(res2,2),unsafe_load(res2,3))
@show resTup2
resfilterd2=filter((x) -> x >0.0 , resTup2)
display(resfilterd2) =#

