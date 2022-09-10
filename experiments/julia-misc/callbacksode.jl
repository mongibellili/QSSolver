using BenchmarkTools
 struct  SSetting
  callbacks :: Vector{Function}
end
function SSetting()
  cbacks=Vector{Function}()
  return SSetting(cbacks)
end

function append_callback(func::Function, ss::SSetting, args::Any...) #:: Function
    cb = ()->func(ss,args...)
    push!(ss.callbacks, cb)
    cb
 end
  
macro callback(expr::Expr)
    expr.head !== :call && error("Expression is not a function call!")
    esc(:(append_callback($(expr.args...))))
end


function integrate( ss::SSetting)
    for callback in ss.callbacks
        callback()
     end
end
function test()
  u=[1.2,2.0]
  du=[0.0,0.0]
  ss=SSetting()
   function odeFunc(ss,du,u)
    du[1]=3*u[2]
    du[2]=u[1]-u[2]
  end
  @callback odeFunc(ss,du,u)
  #integrate(ss)
 # @show du
end
@btime test()