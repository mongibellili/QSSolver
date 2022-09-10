
module m1
        function dothis()
            println("i am doing it")
        end
        macro test2()
            quote
                $(esc(Meta.parse("m1.dothis")))()
            end
        end
end
function test()


 m1.@test2()
 #m1.dothis()
end
test()
