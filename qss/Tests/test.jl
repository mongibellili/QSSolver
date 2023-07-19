function eliminateRef()
    if refEx isa Expr 
      if refEx.args[1]==:+
        a=Symbol((a.args[1]),(refEx.args[2]), "plus",(refEx.args[3]))
      elseif refEx.args[1]==:-
        a=Symbol((a.args[1]),(refEx.args[2]), "minus",(refEx.args[3]))
      end
    else
      a=Symbol((a.args[1]),(refEx))
    end
  end