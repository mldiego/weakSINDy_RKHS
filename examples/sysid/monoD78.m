function dx = monoD(t,x) 
dx(1,1) = 9.9999986920138326240703463554382*x(2) - 10.000002190454324590973556041718*x(1); 
dx(2,1) = 27.999999924693838693201541900635*x(1) - 1.0000036283479403209639713168144*x(2) - 1.0000001146263457485474646091461*x(1)*x(3); 
dx(3,1) = 1.000001022036030917661264538765*x(1)*x(2) - 2.6666629363685387943405658006668*x(3); 
end