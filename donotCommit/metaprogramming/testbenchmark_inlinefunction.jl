
using Markdown
using InteractiveUtils

using PlutoUI

using MacroTools: flatten

using BenchmarkTools
using PlutoUI
using Random
using Polynomials
using PolynomialRoots
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

#with_terminal() do
function test()
	for _ in 1:100
		a = vec(2(rand(1,3) .- 0.5))
       # a = [0.1, -1.0, 2.0]
		#pol = Polynomials.roots(Polynomial(a))
	#p = filter(r->isreal(r) && real(r) > -eps(), pol)
		quadraticsmallestpositiveroot(a[end:-1:1]...)
		#if length(p) !== 0 && n !== nothing 
			#if abs(minimum(real.(p)) - n) > 1.0e-7
			#	@show a pol p n
			#	break
			#end
		#end
	end
end
@btime test()