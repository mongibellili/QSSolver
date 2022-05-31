using TaylorSeries
using Polynomials
using Roots
using NLsolve
using Plots
#using Plots;gr()
using LinearAlgebra
using BenchmarkTools
tₚ = Float64[]
xa₁ = Float64[]
xa₂ = Float64[]
x₁ = Float64[]
x₂ = Float64[]
q₁ = Float64[]
q₂ = Float64[]

function f₁(q::Vector)
  0.01*q[2]
end

function f₂(q::Vector)
  2020.0-100.0*q[1]-100.0*q[2]
end
integrate=TaylorSeries.integrate
function nth_derivative(q₀::Vector{Float64}, res::Vector{Float64}, f::Function, q::Vector{Taylor1{Float64}}, i::Int)
  q[i].coeffs[1:end-1] = q₀

  x = integrate(f(q))
  res[1:end-1] = x.coeffs[2:end-1] - q₀[2:end]
  res[end] = x.coeffs[end]
end

function test(order::Int, Δq::Float64, duration::Float64)
  f = [f₁, f₂]
  q = [Taylor1(zeros(order+1))+0.0, Taylor1(zeros(order+1))+20.0]
  x = [integrate(f[1](q), 0.0), integrate(f[2](q), 20.0)]
  #println(x)
  for i in 1:order-1
    x = [integrate(f[1](x), 0.0), integrate(f[2](x), 20.0)]
   # println(x)
  end
  q = deepcopy(x)
  q[1].coeffs[end] = 0.0
  q[2].coeffs[end] = 0.0
  q[1].coeffs[1] = 0.0 + 0.0*sign(x[1].coeffs[1])*Δq
  q[2].coeffs[1] = 20.0 + 0.0*sign(x[2].coeffs[1])*Δq
  tₓ = [0.0, 0.0]
  tₙ = [0.0, 0.0]
  tₐ = [0.0, 0.0]
  t = 0.0
  tₒ = 0.0
  it = [0, 0]
  start = true
 
  while t < duration && sum(it) < 1000
    t, i = findmin(tₙ)
    it[i] += 1
    
    if t-tₒ > 0.0
      for tt in 0.0:(t-tₒ)/9:t-tₒ
        push!(tₚ, tₒ+tt)
        push!(x₁, evaluate(x[1], tt))#x[1].coeffs[1])
        push!(x₂, evaluate(x[2], tt))#x[2].coeffs[1])
        push!(q₁, evaluate(q[1], tt))#q[1].coeffs[1])
        push!(q₂, evaluate(q[2], tt))#q[2].coeffs[1])

      end
    end
    x[i] = evaluate(x[i], Taylor1([t-tₓ[i], 1.0]))
    tₓ[i] = t
    j = 3 - i
    q[j] = evaluate(q[j], Taylor1([t-tₐ[j], 1.0]))
    #= println("----q_$j: $(q[j])") =#
    tₐ[j] = t
    
    qₙ = deepcopy(q)
    #println(qₙ)
    #println(q===qₙ)
    if x[i].coeffs[end] > 0.0
      qₙ[i] = x[i]+Δq
    else
      qₙ[i] = x[i]-Δq
    end
    qₙ[i].coeffs[end] = 0.0
    xₙ = integrate(f[i](qₙ), x[i].coeffs[1])
#    println("----x_$i: $(x[i])")
#    println("----xₙ: $xₙ")

    if sign(xₙ.coeffs[end]) == sign(x[i].coeffs[end])
      q[i] = deepcopy(qₙ[i])
      #println("if----q_$j: $(q[j])")
      x[i] = deepcopy(xₙ)
      if x[i].coeffs[2:end] == q[i].coeffs[2:end]
        tₙ[i] = Inf
      else
        p = x[i].coeffs-q[i].coeffs
        p[1] = -Δq
        a = fzeros(Polynomial(p),0,5)
        p[1] = Δq
        b = fzeros(Polynomial(p),0,5)
         a = fzeros(Polynomial((x[j]-q[j]-Δq).coeffs),0,5)
         b = fzeros(Polynomial((x[j]-q[j]+Δq).coeffs),0,5)
        tₙ[i] = t + minimum(filter(v->v>0, [a..., b..., Inf]))
        tₙ[j] = t + minimum([filter(v->v>0, fzeros(Polynomial(p),0,5))..., Inf])
      end
    else
      q₀ = deepcopy(q[i].coeffs[1:end-1])
      q₀[1] = deepcopy(x[i].coeffs[1])
      println("else----q2: $(q[2])")
      #q[i].coeffs[1:end-1] = 
      println(nlsolve((qᵢ, res)->nth_derivative(qᵢ, res, f[i], deepcopy(q), i), q₀).zero)
      #println("else----q_$j: $(q[j])")
      println("i= ",i)
      println("after: else----q2: $(q[2])")
      
     # tₙ[i] = Inf
      #x[i] = integrate(f[i](q), xₙ.coeffs[1])
      x[i] = deepcopy(q[i])
      #x[i].coeffs[1] = x̲.coeffs[1]
      x[i].coeffs[1] = xₙ.coeffs[1]
      if abs(q[i].coeffs[1]-x[i].coeffs[1]) > Δq# && tₒ != t
        #q[i].coeffs[1] = x[i].coeffs[1] - sign(q[i].coeffs[1]+x[i].coeffs[1]) * Δq
        #tₙ[i] = t
        x[i].coeffs[1] = q[i].coeffs[1] - 0.5*sign(q[i].coeffs[1]-x[i].coeffs[1]) * Δq
      end
    end
#    println("----q_$i: $(q[i])")
    tₐ[i] = t
    # qₑ = Taylor1[]
    # push!(qₑ, Taylor1([x[1].coeffs..., 0.0]))
    # push!(qₑ, Taylor1([x[2].coeffs..., 0.0]))
    # xₑ = integrate(f[i](qₑ), x[i].coeffs[1])
    # println(xₑ)
    # tₑ = t-xₑ.coeffs[end-1]/xₑ.coeffs[end]/(order+1)
    # println(tₑ)
    # if tₑ > 0 && tₑ < tₙ[i]
    #   tₙ[i] = tₑ
    # end
#     println("----t_$i: $(tₙ[i])")
#     println("----x_$i: $(x[i])")
    for j = 1:2
      x₀ = evaluate(x[j], t-tₓ[j])
      tₓ[j] = t
      #if j==1 println("x[$j].coeffs= ",x[j].coeffs);println(q);println(x₀); end
      x[j] = integrate(f[j](q), x₀)
      #if j==1 println("after: x[$j].coeffs= ",x[j].coeffs) end
#      println("----x_$j: $(x[j])")
      if x[j].coeffs[2:end] == q[j].coeffs[2:end]
        tₙ[j] = Inf
      else
        p = x[j].coeffs-q[j].coeffs
        p[1] = -Δq

        
        a = fzeros(Polynomial(p),0,5)
        p[1] = Δq
        b = fzeros(Polynomial(p),0,5)
        a = fzeros(Polynomial((x[j]-q[j]-Δq).coeffs),0,5) #δq[j]
         b = fzeros(Polynomial((x[j]-q[j]+Δq).coeffs),0,5)
        tₙ[j] = t + minimum(filter(v->v>0, [a..., b..., Inf]))
        tₙ[j] = t + minimum([filter(v->v>0, fzeros(Polynomial(p),0,5))..., Inf])
      end
      #println("----t_$j: $(tₙ[j])")
    end
    tₒ = t
   # println("$t $i $(xₐ[1]) $(abs(x[1].coeffs[1]-xₐ[1])) $(xₐ[2]) $(abs(x[2].coeffs[1]-xₐ[2]))")
  end  #end while
 # println(it)
end
δ = 1e-6
test(4, 0.000005, 5.0)
#= @btime test(4, 0.000005, 5.0) =#
#= plot(tₚ, [abs(x₁-xa₁), abs(x₂-xa₂)]) =#
#display(plot(tₚ, x₁, x₂))
#= display(plot!(tₚ,x₁)) 
display(plot!(tₚ,x₂)) =# 
println("tₚ= ",tₚ)
println("x₁= ",x₁)
println("x₂= ",x₂)
#= println("done") 
readline() =#