function dx = monoD(t,x) 
dx(1,1) = 9.9999999068131728563457727432251*x(2) - 9.9999998352268448797985911369324*x(1); 
dx(2,1) = 27.999999775478499941527843475342*x(1) - 0.99999980878624228353146463632584*x(2) - 0.99999999431531705340603366494179*x(1)*x(3); 
dx(3,1) = 1.0000000308505150314886122941971*x(1)*x(2) - 2.6666667934714496368542313575745*x(3); 
end