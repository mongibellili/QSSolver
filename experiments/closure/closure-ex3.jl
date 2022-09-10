#= function counter()
    x = 0
    () -> x += 1
 end =#
 initialVal=5
 function parentfunction()
    #println(initialVal)
    initial=initialVal
    increm=2
    temp=0
    function counter()
        initial += increm  # initialVal can not be accessed here: a function can only look up once:only from the scope in which the function was defined
        temp=1+initial+temp
    end
    println(temp) # temp here is always 0
    return counter # return a closure: the var initial used inside counter does not get deleted by the GC
 end


i2 = parentfunction() # i2 is the counter function
println( i2()) 

#println(parentfunction())
println( i2())
println( i2())
println( i2())
println( i2())
i3=parentfunction()
println( i3())
println( i3())