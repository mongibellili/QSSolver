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
	@inbounds w = coeffs[1]
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
	@inbounds cohigh = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		colow, cohigh = if colow > eps(F)
            if low > -eps(F)
			    muladd(low, colow, coeff), muladd(high, cohigh, coeff)
            elseif high < eps(F)
                muladd(low, cohigh, coeff), muladd(high, colow, coeff)
            else
                muladd(low, cohigh, coeff), muladd(high, cohigh, coeff)
            end
		elseif cohigh < -eps(F)
            if low > -eps(F)
			    muladd(high, colow, coeff), muladd(low, cohigh, coeff)
            elseif high < eps(F)
                muladd(high, cohigh, coeff), muladd(low, colow, coeff)
            else
                muladd(high, colow, coeff), muladd(low, colow, coeff)
            end
		else
            if low > -eps(F)
			    muladd(high, colow, coeff), muladd(high, cohigh, coeff)
            elseif high < eps(F)
                muladd(low, cohigh, coeff), muladd(low, colow, coeff)
            else
                min(muladd(low, cohigh, coeff), muladd(high, colow, coeff)), 
                max(muladd(low, colow, coeff), muladd(high, cohigh, coeff))
            end
        end
	end
    colow, cohigh
end

@inline function posintervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	@inbounds cohigh = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		colow, cohigh = if colow > eps(F)
			muladd(low, colow, coeff), muladd(high, cohigh, coeff)
		elseif cohigh < -eps(F)
			muladd(high, colow, coeff), muladd(low, cohigh, coeff)
		else
			muladd(high, colow, coeff), muladd(high, cohigh, coeff)
		end
	end
	colow, cohigh
end

@inline function smallestpositiverootintervalnewtonregulafalsi(coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
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
	list = MVector{2N+1, NTuple{2,F}}(undef)
	index = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom′low, codom′high = posintervalhorner(domlow, domhigh, poly′)
		#@show comid codom′low codom′high
		if codom′low < zero(F) < codom′high
			leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
				domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
			elseif comid < -eps(F)
				domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
			else
				return mid
			end
			#@show leftlow lefthigh rightlow righthigh
			if leftlow < lefthigh
				if rightlow < righthigh
					index += 1
					@inbounds list[index] = rightlow, righthigh
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
			if codomlow * codomhigh < eps(F)
				while true
					comid, comid′ = hornerd(mid, poly)
					delta = comid / comid′
					newmid = mid - delta
					if abs(delta) < 1.0e-2mid 
						return newmid
					elseif domlow < newmid < domhigh
						mid = newmid
					else
						if comid * codomlow > -eps(F)
							domlow, codomlow = mid, comid
						elseif comid * codomhigh > -eps(F)
							domhigh, codomhigh = mid, comid
						else
							return mid
						end
						mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
					end
				end
			end
		end
		if index == 0 break end
		@inbounds domlow, domhigh = list[index]
		index -= 1
	end
	return typemax(F)
end

#= using StaticArrays
c=-0.040323840000000166
#b=7597.313022054603-11511.429846502817
b=-3914.116824448214
a=1.021243315770821e8
@show smallestpositiverootintervalnewtonregulafalsi((a, b, c, 1e-6))  =#
#@show smallestpositiverootintervalnewtonregulafalsi((1e-6, c, b, a),pp)
#@show smallestpositiverootintervalnewtonregulafalsi((0.0, -0.2906652780025, -8938.3815305787700, -5.144241008311048e7),pp)