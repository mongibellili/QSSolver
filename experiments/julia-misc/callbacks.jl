using BenchmarkTools
mutable struct  AbstractEvent
    name::String
    starttime::Int
    stoptime::Int
    speed::Float64
    callbacks :: Vector{Function}
end
function append_callback(func::Function, ev::AbstractEvent, args::Any...) #:: Function
    ev.stoptime < 0 && throw(EventProcessed(ev))
    cb = ()->func(ev, args...)
    push!(ev.callbacks, cb)
    cb
  end
  
  macro callback(expr::Expr)
    expr.head !== :call && error("Expression is not a function call!")
    esc(:(append_callback($(expr.args...))))
  end
function computespeed(ev::AbstractEvent,speed::Int)
    #println("stop sim at ", ev.stoptime)
    ev.speed=(ev.stoptime-ev.starttime)/speed
  end
function baserun(ev::AbstractEvent)
    speed=2
    @callback computespeed(ev,speed) 
end
function step(ev::AbstractEvent)
    for callback in ev.callbacks
        callback()
      end
end
function test()
evn1callbacks=Vector{Function}()
event1=AbstractEvent("ev1",1,5,0.0,evn1callbacks)
 baserun(event1)
# display(event1.callbacks[1]);println()
  step(event1)
 # @show event1.speed
end
#@btime test()