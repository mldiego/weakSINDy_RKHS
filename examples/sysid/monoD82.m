function dx = monoD(t,x) 
dx(1,1) = 9.9999984131954988697543740272522*x(2) - 10.000001667969627305865287780762*x(1); 
dx(2,1) = 27.999999635430867783725261688232*x(1) - 1.0000033889525639096973463892937*x(2) - 1.0000000857760369399329647421837*x(1)*x(3); 
dx(3,1) = 1.0000009124578355113044381141663*x(1)*x(2) - 2.6666630451572927995584905147552*x(3); 
end