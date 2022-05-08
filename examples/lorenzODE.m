function dx = lorenzODE(t,x)
dx(1,1) = 10*x(2) - 10*x(1); 
dx(2,1) = 28*x(1) - x(2) -x(1)*x(3); 
dx(3,1) = x(1)*x(2) - 8/3*x(3);
end

