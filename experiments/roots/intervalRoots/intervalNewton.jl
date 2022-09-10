using MacroTools: flatten
using BenchmarkTools

using StaticArrays
@inline function quadraticsmallestpositiveroot(a::F, b::F, c::F) where F <: AbstractFloat
	_a = one(F) / a
	b, c = -0.5b * _a, c * _a
	Δ = muladd(b, b, -c) # b * b - c
	if Δ < -4eps(F) # Complex roots
		typemax(F) 
	elseif Δ > 4eps(F) # Real roots
		if b > eps(F)
			c > eps(F) ? c / (b + sqrt(Δ)) : b + sqrt(Δ)
		elseif b < -eps(F)
			c < eps(F) ? c / (b - sqrt(Δ)) : typemax(F)
		else
			sqrt(-c)
		end
	else # Double real root
		b > -eps(F) ? b : typemax(F) 
	end
 end
@inline function cubicsmallestpositiveroot(a::F, b::F, c::F, d::F) where F <: AbstractFloat
    _a = one(F) / a
    b, c, d = b * _a, c * _a, d * _a
    m = b < c ? b : c
    m = d < m ? d : m
    m > eps(F) && return typemax(F) # Cauchy bound
	_3 = one(F) / 3
	_9 = one(F) / 9
	SQ3 = sqrt(3one(F))
    xₙ = -b * _3
    b²_9 = b * b * _9
    yₙ = muladd(muladd(-2, b²_9, c), xₙ, d)
    δ² = muladd(-_3, c, b²_9)
    h² = 4δ² * δ² * δ²
    Δ = muladd(yₙ, yₙ, -h²)
    if Δ > 4eps(F) # one real root and two complex roots
		p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
		q = δ² / p
		z = xₙ + p + q
		z > -eps(F) ? z : typemax(F)
    elseif Δ < -4eps(F) # three real roots
		θ = abs(yₙ) < eps(F) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
		δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
		z₁ = 2δ * cos(θ)
		z₂ = muladd(-0.5, z₁, xₙ)
		z₃ = SQ3 * δ * sin(θ)
		x₁ = xₙ + z₁
		x₂ = z₂ + z₃
		x₃ = z₂ - z₃
		x = x₁ > -eps(F) ? x₁ : typemax(F)
		x = x₂ > -eps(F) && x₂ < x ? x₂ : x
		x₃ > -eps(F) && x₃ < x ? x₃ : x
    else # double or triple real roots
		δ = cbrt(0.5yₙ)
		x₁ = xₙ + δ
		x₂ = xₙ - 2δ
		x = x₁ > -eps(F) ? x₁ : typemax(F)
		x₂ > -eps(F) && x₂ < x ? x₂ : x
    end
 end


begin
	@inline function quadraticrealroots(b::F, c::F) where F <: AbstractFloat
		b, c = -0.5b, c
		Δ = muladd(b, b, -c) #b * b - c
		if Δ < -4eps(F) # Complex roots
			tuple()
		elseif Δ > 4eps(F) # Real roots
			q = b > 0 ? b + sqrt(Δ) : b - sqrt(Δ)
			(q, c / q)
		else # Double real root
			tuple(b)
		end
	end
	@inline function cubicmaxroot(b::F, c::F, d::F) where F <: AbstractFloat
		_3 = one(F) / 3
		_9 = one(F) / 9
		SQ3 = sqrt(3one(F))
		xₙ = -b * _3
		b²_9 = b * b * _9
		yₙ = muladd(muladd(-2, b²_9, c), xₙ, d) #d + xₙ * (c - 2b²_9)
		δ² = muladd(-_3, c, b²_9) #b²_9 - c *_3
		h² = 4δ² * δ² * δ²
		Δ = muladd(yₙ, yₙ, -h²) #yₙ * yₙ - h²
		if Δ > 4eps(F) # one real root and two complex roots
			p = yₙ < 0 ? cbrt(0.5 * (-yₙ + √Δ)) : cbrt(0.5 * (-yₙ - √Δ))
			q = δ² / p # cbrt(0.5 * (-yₙ - √Δ)) : cbrt(0.5 * (-yₙ + √Δ))
			xₙ + p + q
		elseif Δ < -4eps(F) # three real roots
			θ = abs(yₙ) < eps(F) ? 0.5π * _3 : atan(√abs(Δ) / abs(yₙ)) * _3 # acos(-yₙ / √h²)
			δ = yₙ < 0 ? √abs(δ²) : -√abs(δ²)
			z₁ = 2δ * cos(θ)
			z₂ = muladd(-0.5, z₁, xₙ) #xₙ - 0.5z₁   # x₂ = xₙ + 2δ * cos(2π / 3 - θ)
			z₃ = SQ3 * δ * sin(θ) # x₃ = xₙ + 2δ * cos(4π / 3 - θ)
			x₁ = xₙ + z₁
			x₂ = z₂ + z₃
			x₃ = z₂ - z₃
			x = x₁ > x₂ ? x₁ : x₂
			x₃ > x ? x₃ : x
		else # double or triple real roots
			δ = cbrt(0.5yₙ)
			x₁ = xₙ + δ
			x₂ = xₙ - 2δ
			x₁ > x₂ ? x₁ : x₂
		end
	end
	@inline function quarticsmallestpositiveroot(a::F, b::F, c::F, d::F, e::F) where F <: AbstractFloat
		_a = inv(a)
		b, c, d, e = b * _a, c * _a, d * _a, e * _a
		m = b < c ? b : c
		m = d < m ? d : m
		m = e < m ? e : m
		m > eps(F) && return typemax(F) # Cauchy bound
		xₙ = -0.25b
		b² = b * b
		p = muladd(-0.375, b², c)
		q = muladd(muladd(-0.5, b², 2c), xₙ, d)
		r = muladd(muladd(muladd(-0.1875, b², c), xₙ, d), xₙ, e)
		if abs(r) < 4eps(F)
			xₙ
		elseif abs(q) < 4eps(F)
			y² = quadraticrealroots(p, r)
			l = length(y²)
			l === 0 && return typemax(F)
			y₁ = y²[1] < -eps(F) ? typemax(F) : √y²[1]
			z₁ = xₙ + y₁
			z₂ = xₙ - y₁
			x = z₁ > -eps(F) ? z₁ : typemax(F)
			x = z₂ > -eps(F) && z₂ < x ? z₂ : x
			l === 1 && return x
			y₂ = y²[2] < -eps(F) ? typemax(F) : √y²[2]
			z₃ = xₙ + y₂
			z₄ = xₙ - y₂
			x = z₃ > -eps(F) && z₃ < x ? z₃ : x
			z₄ > -eps(F) && z₄ < x ? z₄ : x
		else
			h² = cubicmaxroot(2p, muladd(p, p, -4r), -q * q)
			h = √h²
			j = 0.5(p + h² - q / h)
			y = quadraticrealroots(h, j)
			l = length(y)
			x = if l === 0
				typemax(F)
			elseif l === 1
				z₁ = xₙ + y[1]
				z₁ > -eps(F) ? z₁ : typemax(F)
			else
				z₁ = xₙ + y[1]
				z₁ = z₁ > -eps(F) ? z₁ : typemax(F)
				z₂ = xₙ + y[2]
				z₂ > -eps(F) && z₂ < z₁ ? z₂ : z₁
			end
			y = quadraticrealroots(-h, r / j)
			l = length(y)
			if l === 0
				x
			elseif l === 1
				z₁ = xₙ + y[1]
				#z₁ > -eps(F) ? z₁ : x
				z₁ > -eps(F) && z₁ < x ? z₁ : x
			else
				z₁ = xₙ + y[1]
				#x = z₁ > -eps(F) ? z₁ : x
				x = z₁ > -eps(F) && z₁ < x ? z₁ : x
				z₂ = xₙ + y[2]
				z₂ > -eps(F) && z₂ < x ? z₂ : x
			end
		end
	end
end
@generated function smallest(arg1::F, arg2::F, args::F...) where F <: AbstractFloat
	ex = quote m = arg1 < arg2 ? arg1 : arg2 end
	for i in 1:length(args)
		ex = quote
			$ex
			@inbounds m = args[$i] < m ? args[$i] : m
		end
	end
	ex |> flatten
end
@generated function horner(x::Float64, coeff1::Float64, coeffs::Float64...)::Float64# where F <: AbstractFloat
 	l = length(coeffs)
	l === 0 && return quote coeff1 end
	ex = :(coeff1)
	for i in 1:l
		ex = :(@inbounds muladd(x, $ex, coeffs[$i]))
	end
	ex |> flatten
end
@generated function posintervalhorner(low::Float64, high::Float64, coeff1::Float64, coeffs::Float64...)::Float64# where F <: AbstractFloat
	l = length(coeffs)
	l === 0 && return quote coeff1, coeff1 end
	ex = quote colow, cohigh = if coeff1 > eps()
				@inbounds muladd(low, coeff1, coeffs[1]), muladd(high, coeff1, coeffs[1])
			elseif coeff1 < -eps()
				@inbounds muladd(high, coeff1, coeffs[1]), muladd(low, coeff1, coeffs[1])
			else
				@inbounds muladd(high, coeff1, coeffs[1]), muladd(high, coeff1, coeffs[1])
			end
		 end
	for i in 2:l
	ex = quote
		$ex
		colow, cohigh = if colow > eps()
				@inbounds muladd(low, colow, coeffs[$i]), muladd(high, cohigh, coeffs[$i])
			elseif cohigh < -eps()
				@inbounds muladd(high, colow, coeffs[$i]), muladd(low, cohigh, coeffs[$i])
			else
				@inbounds muladd(high, colow, coeffs[$i]), muladd(high, cohigh, coeffs[$i])
			end
		end
	end
	ex |> flatten
end
struct List{F<:AbstractFloat}
	prev :: Union{Nothing, List{F}}
	low :: F
	high :: F
end
#= domlow=0.0
domhigh=1.23
dercoeff1=3.0
dercoeff2=0.78
dercoeff3=2.2
coeff4=1.056
@btime horner(domhigh, dercoeff1,dercoeff2,dercoeff3,coeff4)#5.759 ns (0 allocations: 0 bytes)
@btime posintervalhorner(domlow, domhigh, dercoeff1,dercoeff2,dercoeff3)# 7.351 ns (0 allocations: 0 bytes) =#
 
@generated function gensmallestpositiverootintervalnewton(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		mm > eps(F) && return typemax(F)
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true    # instead of starting with 1 point, find the next hope (root) and iterate to a better hope, we start with 2pts (interval) and try to tighten the interval hopefully contain the root
			@show domlow domhigh
			mid = 0.5(domlow + domhigh)
            @show mid
			comid = horner(mid, $(syms...))  # f(middle_of_interval)
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))  #calculate mixed derivatives(poshornerinterval) ...of the interval bounds
			@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)  # if mixed der of interval bounds are of opposite sign ie large difference of der (except one case) ie bounds are away from each other
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)  #find/calculate new better bounds for the interval
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid # this is the special case where mixed der of different signs can yield a sol (mid is the root)
				end
                @show leftlow, lefthigh, rightlow, righthigh
				if !(domhigh < leftlow || lefthigh < domlow)  # if lh (ie new lower bound) is better (>) than old leftbound (domlow) # we want to update the .....surprisingly the upperbound!!!!
					if !(domhigh < rightlow || righthigh < domlow) # same thing if rl (new higher bound) is better (<) than old rightbound (domhigh)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)  #we store only the lowerbound update in a list if both better bounds
                       @show list
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh # lowerbound stayed the same(leftlow is -Inf) and upperbound updated with the better one. why only rightbound???
					println("update domlow,domhigh= ",domlow,domhigh)
				elseif !(domhigh < rightlow || righthigh < domlow)# same thing we check right bound: if a better (<) higher bound, we want to update 
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh# here we update the lowerbound, the upper bound is kept the same (righthigh is Inf)
                    println("new domlow domhigh = ", domlow, domhigh)
				elseif list === nothing  # if both new bounds are worse than old bounds and no previous info
					return typemax(F)
				else # both new bounds are worse than old bounds and there is previous info
					list, domlow, domhigh = list.prev, list.low, list.high # go back to previous state
                    println("update from list: list gets prev, domlow,domhigh= ",domlow,domhigh)
				end
			else  #mixed der of the bounds are of same sign ie they are getting close and means the bounds themselves are getting close
				low, high = if comid > eps(F) # calculate new bounds
					mid - comid / codom′low, mid - comid  / codom′high
				elseif comid < -eps(F)
					mid - comid  / codom′high, mid - comid / codom′low
				else
					return mid
				end
                println(" if codomlow and codomhigh are of same sign; low,high = ",low,high)
				low = domlow < low ? low : domlow #update lowerbound if better
				high = domhigh > high ? high : domhigh #update higher bound if better
				@show low high
				domlow, domhigh = if high - low > -eps(F)  # we still have a valid interval in which the possible root might be found
					if high - low < 1e-8mid  # if interval tight enough, root is found, return 
                        println("very tight interval")
						return 0.5(low + high)
					end
					low, high
				elseif list === nothing # we do not have a valid interval and no previous info saved
					return typemax(F)
				else # not valid interval and we have previous info
					list, low, high = list.prev, list.low, list.high #go back to previous state
					low, high
				end
                println(" if codomlow and codomhigh are of same sign; domlow,domhigh = ",domlow,domhigh)
			end
		end
	end
	ex |> flatten
end

#= function gensmallestpositiverootintervalnewton(coeff1::Float64,coeff2::Float64,coeff3::Float64,coeff4::Float64 )::Float64# where F <: AbstractFloat
	
	_coeff1 = 1.0/coeff1
	coeff1 = 1.0
	coeff2 = _coeff1 * coeff2
	coeff3 = _coeff1 * coeff3
	coeff4 = _coeff1 * coeff4
	dercoeff1 = 3.0
	dercoeff2 = 2coeff2
	dercoeff3 = coeff3
		mm = smallest(coeff2,coeff3,coeff4)
		#mm > eps() && return Inf
		mm > 0.0 && return Inf
		domlow, domhigh = 0.0, 1.0 - mm
		list = nothing
		while true    # instead of starting with 1 point, find the next hope (root) and iterate to a better hope, we start with 2pts (interval) and try to tighten the interval hopefully contain the root
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
            #@show mid
			comid = horner(mid, coeff1,coeff2,coeff3,coeff4)  # f(middle_of_interval)
			codom′low, codom′high = posintervalhorner(domlow, domhigh, dercoeff1,dercoeff2,dercoeff3)  #calculate mixed derivatives(poshornerinterval) ...of the interval bounds
			#@show comid codom′low codom′high
			if codom′low < -eps() && codom′high > eps()  # if mixed der of interval bounds are of opposite sign ie large difference of der (except one case) ie bounds are away from each other
				leftlow, lefthigh, rightlow, righthigh = if comid > eps()  #find/calculate new better bounds for the interval
					-Inf, mid - comid / codom′high, mid - comid / codom′low, Inf
				elseif comid < -eps()
					-Inf, mid - comid / codom′low, mid - comid / codom′high, Inf
				else
					return mid # this is the special case where mixed der of different signs can yield a sol (mid is the root)
				end
                #@show leftlow, lefthigh, rightlow, righthigh
				if !(domhigh < leftlow || lefthigh < domlow)  # if lh (ie new lower bound) is better (>) than old leftbound (domlow) # we want to update the .....surprisingly the upperbound!!!!
					if !(domhigh < rightlow || righthigh < domlow) # same thing if rl (new higher bound) is better (<) than old rightbound (domhigh)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)  #we store only the lowerbound update in a list if both better bounds
                       #@show list
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh # lowerbound stayed the same(leftlow is -Inf) and upperbound updated with the better one. why only rightbound???
					#println("update domlow,domhigh= ",domlow,domhigh)
				elseif !(domhigh < rightlow || righthigh < domlow)# same thing we check right bound: if a better (<) higher bound, we want to update 
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh# here we update the lowerbound, the upper bound is kept the same (righthigh is Inf)
                   # println("new domlow domhigh = ", domlow, domhigh)
				elseif list === nothing  # if both new bounds are worse than old bounds and no previous info
					return Inf
				else # both new bounds are worse than old bounds and there is previous info
					list, domlow, domhigh = list.prev, list.low, list.high # go back to previous state
                    #println("update from list: list gets prev, domlow,domhigh= ",domlow,domhigh)
				end
			else  #mixed der of the bounds are of same sign ie they are getting close and means the bounds themselves are getting close
				low, high = if comid > eps() # calculate new bounds
					mid - comid / codom′low, mid - comid  / codom′high
				elseif comid < -eps()
					mid - comid  / codom′high, mid - comid / codom′low
				else
					return mid
				end
                #println(" if codomlow and codomhigh are of same sign; low,high = ",low,high)
				low = domlow < low ? low : domlow #update lowerbound if better
				high = domhigh > high ? high : domhigh #update higher bound if better
				#@show low high
				domlow, domhigh = if high - low > -eps()  # we still have a valid interval in which the possible root might be found
					if high - low < 1e-8mid  # if interval tight enough, root is found, return 
                        #println("very tight interval")
						return 0.5(low + high)
					end
					low, high
				elseif list === nothing # we do not have a valid interval and no previous info saved
					return Inf
				else # not valid interval and we have previous info
					list, low, high = list.prev, list.low, list.high #go back to previous state
					low, high
				end
                #println(" if codomlow and codomhigh are of same sign; domlow,domhigh = ",domlow,domhigh)
			end
		end
  return -1
	
end

#@btime gensmallestpositiverootintervalnewton(1.0, -2.0, -1.0)#64.912 ns (1 allocation: 32 bytes)
#@show gensmallestpositiverootintervalnewton(1.0, -2.0, -1.0)#2.414213562506035#@show quadraticsmallestpositiveroot(1.0, -2.0, -1.0)#1.649 ns (0 allocations: 0 bytes)
#@btime gensmallestpositiverootintervalnewton(-1.33,2.0,1.0,6.0)#126.486 ns (2 allocations: 64 bytes)....2.5155562382541037...cubic:18ns
#@btime gensmallestpositiverootintervalnewton(1.33,2.0,-2.63,6.0)#32.639 ns (0 allocations: 0 bytes) ....Inf
#@show gensmallestpositiverootintervalnewton(2.0,-3.0,-0.5,-5.0)
#@btime gensmallestpositiverootintervalnewton(-1.33,2.0,1.0,6.0,2.2)#171.461 ns (2 allocations: 64 bytes)
#@btime gensmallestpositiverootintervalnewton(-1.33,2.0,1.0,6.0,2.2,3.14)#220.877 ns (2 allocations: 64 bytes)
@generated function smallestpositiverootintervalnewtonpure(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F) # now instead of highbound>lowbound (ie valid interval), we here actually check if the interval contains the root (we use sym not sym')
					x = mid
					f = comid
					while true   # this while loop is pure newton
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						if abs(delta) < 1.0e-8x return x - delta end
						x -= delta
						f = horner(x, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end
#@btime smallestpositiverootintervalnewtonpure(1.0, -2.0, -1.0)#50.002 ns (1 allocation: 32 bytes)# 
#@btime smallestpositiverootintervalnewtonpure(-1.33,2.0,1.0,6.0)# 94.988 ns (2 allocations: 64 bytes)
#@btime smallestpositiverootintervalnewtonpure(1.33,2.0,-2.63,6.0)#32.019 ns (0 allocations: 0 bytes)
#@btime smallestpositiverootintervalnewtonpure(-1.33,2.0,1.0,6.0,2.2)#118.146 ns (2 allocations: 64 bytes)...quartic
#@btime smallestpositiverootintervalnewtonpure(-1.33,2.0,1.0,6.0,2.2,3.14)#166.758 ns (2 allocations: 64 bytes) ...quintic...about 74µs for recompute for ft=10
@generated function smallestpositiverootintervalnewtonhalley(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	syms′′ = [gensym() for _ in 1:l-1]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	for i in 1:l-1
		ex = quote
			$ex
			$(syms′′[i]) = $(l-i) * $(syms′[i])
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						f′′ = horner(x, $(syms′′...))
						delta = f / f′
						delta = delta / (one(F) - 0.5delta * f′′ / f′)
						if abs(delta) < 1.0e-8x return x - delta end
						x -= delta
						f = horner(x, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end
#@btime smallestpositiverootintervalnewtonhalley(1.0, -2.0, -1.0)#82.505 ns (1 allocation: 32 bytes)
#@btime smallestpositiverootintervalnewtonhalley(-1.33,2.0,1.0,6.0)#   139.965 ns (2 allocations: 64 bytes)
#@btime smallestpositiverootintervalnewtonhalley(-1.33,2.0,1.0,6.0,2.2)#155.568 ns (2 allocations: 64 bytes)
#@btime smallestpositiverootintervalnewtonhalley(-1.33,2.0,1.0,6.0,2.2,3.14)#195.203 ns (2 allocations: 64 bytes)
@generated function smallestpositiverootintervalnewtonregulafalsi(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				if codomlow * codomhigh < eps(F)
					side = 0
					while domhigh - domlow > 0.5e-8(domhigh + domlow)# now instead of leaving the loop on the condition xn close to xn+1, we leave when interval is small
						#@show domlow domhigh codomlow codomhigh
						mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow)# usually mid is 0.5(domlow + domhigh), here we find the secant (equation of a line=0)
						comid = horner(mid, $(syms...))
						if comid * codomlow < -eps(F) # which point to pick ; if colow and comid are opposite, make comid the new cohi
							domhigh = mid
							codomhigh = comid
							if side === -1   # if the secant still on the left, 
								codomlow *= 0.5#make colow smaller so that the secant gets closer to the root faster, avoids regularfalsi failures
							end
							side = -1
						elseif comid * codomhigh < -eps(F)
							domlow = mid
							codomlow = comid
							if side === 1
								codomhigh *= 0.5
							end
							side = 1
						else
							break
						end
					end
					return 0.5(domhigh + domlow)
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end
#@show smallestpositiverootintervalnewtonregulafalsi(-1.33,2.0,1.0,6.0)# 2.5155562373313574
#@btime smallestpositiverootintervalnewtonregulafalsi(-1.33,2.0,1.0,6.0)# 2.5155562373313574# 142.592 ns (2 allocations: 64 bytes)

@generated function smallestpositiverootintervalnewtonridders(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				if codomlow * codomhigh < eps(F)
					ans = typemax(F)
					while true
						s = sqrt(comid*comid - codomlow*codomhigh)
						#if s === 0.0 return ans end
						next = mid + (codomlow > codomhigh ? 1 : -1) * (mid-domlow) * comid / s
						if abs(ans - next) < 1e-7ans return ans end
						ans = next
						conext = horner(next, $(syms...))
						if comid * conext < 0
							domlow = mid
							codomlow = comid
							domhigh = next
							codomhigh = conext
						elseif codomlow * conext < 0
							domhigh = next
							codomhigh = conext
						else
							domlow = next
							codomlow = conext
						end
						if abs(domhigh - domlow) < 1e-8ans return ans end
						mid = 0.5(domlow + domhigh)
						comid = horner(mid, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end =#

@generated function smallestpositiverootintervalnewtonrobust(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f*codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f*codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							x = 0.5(domlow + domhigh)
						end
						f = horner(x, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end


@generated function smallestpositiverootintervalnewtonrobustregulafalsi(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = nothing
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   list = List(list, domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f*codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f*codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							newx = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							x = if domlow < newx < domhigh   # there is a bad case where it won't get here ...domlow and domhigh to the right of the root and the curve is almost flat to the right of the root. 
								newx
							else
								0.5(domlow + domhigh) # bisection
							end
						end
						f = horner(x, $(syms...))
					end
				elseif list === nothing
					return typemax(F)
				else
					list, domlow, domhigh = list.prev, list.low, list.high
				end
			end
		end
	end
	ex |> flatten
end
#@show smallestpositiverootintervalnewtonrobustregulafalsi(1.0, -4.0, 3.9999999)#

# ╔═╡ 46820ceb-a88b-4494-9bc3-2de11ecd49e7
#= @generated function smallestpositiverootintervalnewtonrobustregulafalsistatic(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = MVector{$(2l), NTuple{2, F}}(undef)
		index = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
						index += 1
					   	@inbounds list[index] = tuple(domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index]
					index -= 1
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f * codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f * codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							newx = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							x = domlow < newx < domhigh ? newx : 0.5(domlow + domhigh) # bisection
						end
						f = horner(x, $(syms...))
					end
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index]
					index -= 1
				end
			end
		end
	end
	ex |> flatten
end
@show smallestpositiverootintervalnewtonrobustregulafalsistatic(1.0, -4.0, 3.9)# =#

#= 
@generated function smallestpositiverootintervalnewtonrobustregulafalsiarray(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = MVector{40, F}(undef)
		index = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
						index += 1
					   	@inbounds list[index] = domlow < rightlow ? rightlow : domlow
						index += 1
					   	@inbounds list[index] = domhigh < righthigh ? domhigh : righthigh
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif index === 0
					return typemax(F)
				else
					@inbounds domhigh = list[index]
					index -= 1
					@inbounds domlow = list[index]
					index -= 1
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f*codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f*codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							newx = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							x = if domlow < newx < domhigh
								newx
							else
								0.5(domlow + domhigh) # bisection
							end
						end
						f = horner(x, $(syms...))
					end
				elseif index === 0
					return typemax(F)
				else
					@inbounds domhigh = list[index]
					index -= 1
					@inbounds domlow = list[index]
					index -= 1
				end
			end
		end
	end
	ex |> flatten
end
 =#
#= @generated function gensmallestpositiverootintervalnewtonstatic(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		mm > eps(F) && return typemax(F)
		domlow, domhigh = zero(F), one(F) - mm
		list = MVector{20, NTuple{2, F}}(undef)
		index = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
						index += 1
					   	@inbounds list[index] = tuple(domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index]
					index -= 1
				end
			else
				low, high = if comid > eps(F)
					mid - comid / codom′low, mid - comid  / codom′high
				elseif comid < -eps(F)
					mid - comid  / codom′high, mid - comid / codom′low
				else
					return mid
				end
				low = domlow < low ? low : domlow
				high = domhigh > high ? high : domhigh
				#@show low high
				if high - low > -eps(F)
					if high - low < 1e-8mid
						return 0.5(low + high)
					end
					domlow, domhigh = low, high
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index]
					index -= 1
				end
			end
		end
	end
	ex |> flatten
end =#

@generated function smallestpositiverootintervalnewtonrobustregulafalsistaticref(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = (Ref{NTuple{2,F}}(), Ref{NTuple{2,F}}(), Ref{NTuple{2,F}}(), Ref{NTuple{2,F}}(), Ref{NTuple{2,F}}(), Ref{NTuple{2,F}}())
		index = 0
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
						index += 1
					   	@inbounds list[index][] = tuple(domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh)
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index][]
					index -= 1
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f * codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f * codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							newx = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							x = domlow < newx < domhigh ? newx : 0.5(domlow + domhigh) # bisection
						end
						f = horner(x, $(syms...))
					end
				elseif index === 0
					return typemax(F)
				else
					@inbounds domlow, domhigh = list[index][]
					index -= 1
				end
			end
		end
	end
	ex |> flatten
end


# ╔═╡ 37dc92fc-f7de-4569-b4c3-25ce858add1d
mutable struct Tupple2FloatList{N, F<:AbstractFloat}
	data :: NTuple{N, NTuple{2, F}}
	index :: Int
	function Tupple2FloatList{N, F}() where {N, F<:AbstractFloat}
		t2fl = new{N, F}()
		t2fl.index = 0
		t2fl
	end
end

# ╔═╡ ca5fde85-5b99-4b55-a520-07d0a7e670a7
@inline function Base.push!(t2fl::Tupple2FloatList{N, F}, val::NTuple{2, F}) where {N, F<:AbstractFloat}
	i = t2fl.index + 1
	GC.@preserve t2fl unsafe_store!(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(t2fl)), convert(NTuple{2, F}, val), i)
	t2fl.index = i
	t2fl
end

# ╔═╡ 093a262a-3345-4180-b042-def2120b7fb6
@inline function Base.pop!(t2fl::Tupple2FloatList{N, F} ) where {N, F<:AbstractFloat}
	i = t2fl.index
	val = GC.@preserve t2fl unsafe_load(Base.unsafe_convert(Ptr{NTuple{2, F}}, pointer_from_objref(t2fl)), i)
	t2fl.index = i - 1
	val
end

# ╔═╡ 10f6769e-5386-476d-87eb-ac3ec97d0cf4
@inline function Base.length(t2fl::Tupple2FloatList{N, F} ) where {N, F<:AbstractFloat}
	t2fl.index
end

# ╔═╡ 198031bc-48c4-48aa-8bc2-18ff732e5b2a
@generated function smallestpositiverootintervalnewtonrobustregulafalsit2fl(coeff1::F, coeffs::F...) where F <: AbstractFloat
	l = length(coeffs)
	syms = [gensym() for _ in 1:l+1]
  	syms′ = [gensym() for _ in 1:l]
	ex = quote
		_coeff1 = inv(coeff1)
		$(syms[1]) = one(F)
	end
	for i in 1:l
		ex = quote
			$ex
			$(syms′[i]) = $(l-i+1) * $(syms[i])
			@inbounds $(syms[i+1]) = _coeff1 * coeffs[$i]
		end
	end
	ex = quote
		$ex
		mm = smallest($(syms[2:l+1]...))
		if mm > eps(F) return typemax(F) end
		domlow, domhigh = zero(F), one(F) - mm
		list = Tupple2FloatList{$(2l), F}()
		while true
			#@show domlow domhigh
			mid = 0.5(domlow + domhigh)
			comid = horner(mid, $(syms...))
			codom′low, codom′high = posintervalhorner(domlow, domhigh, $(syms′...))
			#@show comid codom′low codom′high
			if codom′low < -eps(F) && codom′high > eps(F)
				leftlow, lefthigh, rightlow, righthigh = if comid > eps(F)
					typemin(F), mid - comid / codom′high, mid - comid / codom′low, typemax(F)
				elseif comid < -eps(F)
					typemin(F), mid - comid / codom′low, mid - comid / codom′high, typemax(F)
				else
					return mid
				end
				if !(domhigh < leftlow || lefthigh < domlow)
					if !(domhigh < rightlow || righthigh < domlow)
					   	push!(list, (domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh))
					end
					domlow, domhigh = domlow < leftlow ? leftlow : domlow, domhigh < lefthigh ? domhigh : lefthigh
				elseif !(domhigh < rightlow || righthigh < domlow)
					domlow, domhigh = domlow < rightlow ? rightlow : domlow, domhigh < righthigh ? domhigh : righthigh
				elseif length(list) === 0
					return typemax(F)
				else
					domlow, domhigh = pop!(list)
				end
			else
				codomlow = horner(domlow, $(syms...))
				codomhigh = horner(domhigh, $(syms...))
				#@show domlow domhigh codomlow codomhigh
				if codomlow * codomhigh < eps(F)
					x = mid
					f = comid
					while true
						f′ = horner(x, $(syms′...))
						#@show x f f′
						delta = f / f′
						newx = x - delta
						if abs(delta) < 1.0e-8x 
							return newx
						elseif domlow < newx < domhigh
							x = newx
						else
							if f * codomlow > -eps(F)
								domlow = x
								codomlow = f
							elseif f * codomhigh > -eps(F)
								domhigh = x
								codomhigh = f
							else
								return x
							end
							newx = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regulafalsi
							x = domlow < newx < domhigh ? newx : 0.5(domlow + domhigh) # bisection
						end
						f = horner(x, $(syms...))
					end
				elseif length(list) === 0
					return typemax(F)
				else
					domlow, domhigh = pop!(list)
				end
			end
		end
	end
	ex |> flatten
end

@btime smallestpositiverootintervalnewtonrobustregulafalsit2fl(-1.33,2.0,1.0,6.0)# 99.994 ns (0 allocations: 0 bytes)....2.5155562382541037...cubic:18ns