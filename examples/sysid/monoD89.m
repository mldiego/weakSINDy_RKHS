function dx = monoD(t,x) 
dx(1,1) = 9.9999999179799488047137856483459*x(2) - 9.9999998279181454563513398170471*x(1); 
dx(2,1) = 27.999999952458892948925495147705*x(1) - 0.99999984430496624554507434368134*x(2) - 0.99999999989802290656371042132378*x(1)*x(3); 
dx(3,1) = 1.0000000378890945285093039274216*x(1)*x(2) - 2.6666668619741358270402997732162*x(3); 
end