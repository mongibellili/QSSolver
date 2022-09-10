integratorCache=zcf[i](q,d,initTime + t,taylorOpsCache) #test this evaluation
oldsignValue[i,2]=integratorCache.coeffs[1] #value
oldsignValue[i,1]=sign(integratorCache.coeffs[1]) #sign modify 
computeNextEventTime(i,integratorCache,oldsignValue,initTime,  nextEventTime, quantum)