#= function counter()
    x = 0
    () -> x += 1
 end =#
 initialVal=5
 function parentfunction()
    #println(initialVal)
    initial=initialVal

    function counter()
        initial += 1  
       
    end
    
    return counter # return a closure: the var initial used inside counter does not get deleted by the GC
 end


i2 = parentfunction() # i2 is the counter function
println( i2()) 

println(parentfunction())
println( i2())
println( i2())

parentfunction()