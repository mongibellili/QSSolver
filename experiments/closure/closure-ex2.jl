#= 
 function parentfunction(x)
   
    function counter()
        x += 1  
    end
  
    return counter # 
 end

initVal=0
i2 = parentfunction(initVal) # i2 is the counter function
println( i2()) 

#println(parentfunction())
println( i2())
println( i2())
initVal=0
println( i2())
println( i2())
i3=parentfunction(initVal)
println( i3())
println( i3()) =#

function parentfunction()
   x=0
    function counter()
        x += 1  
    end
  
    return counter # 
 end

initVal=0
i2 = parentfunction() # i2 is the counter function
println( i2()) 

#println(parentfunction())
println( i2())
println( i2())
initVal=0
println( i2())
println( i2())
i3=parentfunction()
println( i3())
println( i3())