function dx = monoD(t,x) 
dx(1,1) = 9.9999997754948708461597561836243*x(2) - 9.9999999357914930442348122596741*x(1); 
dx(2,1) = 27.999992959164956118911504745483*x(1) - 0.99999783949124321225099265575409*x(2) - 0.99999972607133713609073311090469*x(1)*x(3); 
dx(3,1) = 1.0000002162355485779698938131332*x(1)*x(2) - 2.6666698004573845537379384040833*x(3); 
end