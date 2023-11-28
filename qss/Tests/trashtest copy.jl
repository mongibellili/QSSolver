
setprecision(BigFloat,90)

aii, aij, aji, ajj, h, xi, xjaux, uij, uij2, uji, uji2 = BigFloat(-1.0e6), BigFloat(1000.0), BigFloat(1.0e6), BigFloat(-1000.0157261488598), BigFloat(8.858170228361118), BigFloat(0.0009738922631926892), BigFloat(0.9538344383858633), BigFloat(0.024830940874835505), BigFloat(-0.02117642588838642), BigFloat(-5.135234459885396e-10), BigFloat(0.00014385727589605324)

#aii, aij, aji, ajj, h, xi, xjaux, uij, uij2, uji, uji2 = -1.0e6, 1000.0, 1.0e6, -1000.0157261488598, 8.858170228361118, 0.0009738922631926892, 0.9538344383858633, 0.024830940874835505, -0.02117642588838642, -5.135234459885396e-10, 0.00014385727589605324
h_2=h*h;h_3=h_2*h;h_4=h_3*h;h_5=h_4*h;h_6=h_5*h
aiijj=aii+ajj
aiijj_=aij*aji-aii*ajj
#Δ1=(1-h*aii)*(1-h*ajj)-h*h*aij*aji
Δ1=1.0-h*(aiijj)-h_2*(aiijj_)
if abs(Δ1)==0.0
  Δ1=1e-30
  @show Δ1
end

Δ1_2=Δ1*Δ1
Δ1_2_=1.0-2.0*h*aiijj+h_2*(aiijj*aiijj-2.0*aiijj_)+2.0*h_3*aiijj*aiijj_+h_4*aiijj_*aiijj_
#= if abs(Δ1_2-Δ1_2_) >1e-3
  @show Δ1_2,Δ1_2_
end =#
αii=(aii*(1.0-h*ajj)+h*aij*aji)/Δ1
# αij=((1-h*ajj)*aij+h*aij*ajj)/Δ1
αij=aij/Δ1
αji=aji/Δ1
αjj=(h*aji*aij+(1.0-h*aii)*ajj)/Δ1
βii_=1.0+h*(αii-aii)-h_2*(aii*αii+aij*αji)/2.0
βii=1.0+(h_2*(aij*aji+aii*aii)/(2.0*Δ1)+h_3*aii*(aiijj_)/(2.0*Δ1))
βij_=h*(αij-aij)-h_2*(aii*αij+aij*αjj)/2.0
βij=h_2*aij*(aiijj)/(2.0*Δ1)+h_3*aij*( aiijj_)/(2.0*Δ1)
βji_=h*(αji-aji)-h_2*(aji*αii+ajj*αji)/2.0
βji=h_2*aji*(aiijj)/(2.0*Δ1)+h_3*aji*(aiijj_)/(2.0*Δ1)

βjj_=1.0+h*(αjj-ajj)-h_2*(aji*αij+ajj*αjj)/2.0
#β__jj=1+(h*h*(aji*aij+ajj*ajj)/2+h*h*h*(ajj*aij*aji-aii*ajj*ajj)/2)/Δ1
βjj=1.0+(h_2*(aji*aij+ajj*ajj)/(2.0*Δ1)+h_3*ajj*( aiijj_)/(2.0*Δ1))


Δ22=βii*βjj-βij*βji


Δ2=1.0+(h_2*(aii*aii+ajj*ajj+2.0*aij*aji)+h_3*(aiijj)*(aiijj_))/(2.0*Δ1)+(h_4*(aiijj_)*(aiijj_)-h_5*(aiijj)*(aiijj_)*(aiijj_)-h_6*aiijj_*aiijj_*aiijj_)/(4.0*Δ1*Δ1)
Δ2_=(4.0*Δ1*Δ1+2.0*Δ1*(h_2*(aii*aii+ajj*ajj+2.0*aij*aji)+h_3*(aiijj)*(aiijj_))+h_4*(aiijj_)*(aiijj_)-h_5*(aiijj)*(aiijj_)*(aiijj_)-h_6*aiijj_*aiijj_*aiijj_)/(4.0*Δ1*Δ1)

#my expres
Δ2__=(4.0-8.0*h*aiijj+h_2*(6.0*aiijj*aiijj-4.0*aiijj_)+h_3*(6.0*aiijj*aiijj_-2.0*aiijj*aiijj*aiijj)+h_4*(aiijj_*aiijj_-4.0*aiijj_*aiijj*aiijj)-3.0*h_5*aiijj*aiijj_*aiijj_-h_6*aiijj_*aiijj_*aiijj_)/(4.0*Δ1*Δ1) 


α=aii*aii+ajj*ajj+2.0*aij*aji
β=aiijj*aiijj_



#correct
ΔΔ=(4.0-8.0*h*aiijj+4.0*h*h*(aiijj*aiijj-2.0*aiijj_)+8.0*h_3*aiijj*aiijj_+4.0*h_4*aiijj_*aiijj_+2.0*α*h*h-2*α*aiijj*h_3-2.0α*h_4*aiijj_+2.0*β*h_3-2.0*h_4*aiijj*β-3.0*h_5*aiijj*aiijj_*aiijj_+h_4*aiijj_*aiijj_-h_6*aiijj_*aiijj_*aiijj_)/(4.0*Δ1*Δ1) 
  @show Δ2,Δ2_,Δ2__,ΔΔ
 

 # Δ2__squar=(h_2*(6.0*aiijj*aiijj-4.0*aiijj_))/(4.0*Δ1*Δ1) 
 # ΔΔ_square=(4.0*h*h*(aiijj*aiijj-2.0*aiijj_)+2.0*α*h*h)/(4.0*Δ1*Δ1) 
 # @show Δ2__squar,ΔΔ_square

#=  ΔΔ_cub=(8.0*h_3*aiijj*aiijj_-2*α*aiijj*h_3+2.0*β*h_3)/(4.0*Δ1*Δ1) 
  Δ2__cub=h_3*(6.0*aiijj*aiijj_-2.0*aiijj*aiijj*aiijj)/(4.0*Δ1*Δ1) 
  @show ΔΔ_cub,Δ2__cub =#


 ΔΔ_quart=(4.0*h_4*aiijj_*aiijj_-2.0α*h_4*aiijj_-2.0*h_4*aiijj*β+h_4*aiijj_*aiijj_)/(4.0*Δ1*Δ1) 
 mydelta_quart=(h_4*(aiijj_*aiijj_-4.0*aiijj_*aiijj*aiijj))/(4.0*Δ1*Δ1) 
 @show mydelta_quart,ΔΔ_quart