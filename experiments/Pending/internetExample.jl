function permutations(array,n)
    exp = quote
        l = []
        @nloops $n i _ -> 1:length($array) begin
            idx =[@ntuple($n,i)...]
            push!(l,$array[idx])
        end
        l
    end
    eval(exp)
end

permutations("ACTG",2)