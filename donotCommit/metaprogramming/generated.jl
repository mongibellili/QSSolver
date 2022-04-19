@generated function horner(x::F, coeff1::F, coeffs::F...) where F <: AbstractFloat
    l = length(coeffs)
   l === 0 && return quote coeff1 end
   ex = :(coeff1)
   for i in 1:l
       ex = :(@inbounds muladd(x, $ex, coeffs[$i]))
   end
   ex |> flatten
end

# ╔═╡ 439e5aea-5a91-11eb-30e9-63c970d2b687
with_terminal() do
   @code_warntype horner(2.0, -1.0, 2.0, 3.0, 4.0)
end