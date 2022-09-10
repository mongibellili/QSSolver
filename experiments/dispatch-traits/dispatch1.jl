displaypercent(x)=println("$(x*100) %")
displaypercent(x::Rational)=println(round(100x; digits=2)," %")
displaypercent(x::String)=println("$x %")


function display_percent(str::AbstractString)
   # if occursin(r"\d+%$", str)  # any combination of numbers, followed by a percent sign
    if occursin(r"^\d*\.?\d+\s*%$", str)
        println(str)
    else
        throw(DomainError(str, "Not valid percentage format"))
    end
end


#display_percent("9.998%")
abstract type  Fraction end
struct Half<:Fraction end
displaypercent(x::Half)=println("50%")
#displaypercent(x::Fraction)=println("50%")
h=Half()
displaypercent(Half())

