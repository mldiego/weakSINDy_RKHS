function dx = monoD(t,x) 
dx(1,1) = 9.9999996684036887018010020256042*x(2) - 9.9999998775983840459957718849182*x(1); 
dx(2,1) = 27.999992561057297280058264732361*x(1) - 0.9999975649185444126487709581852*x(2) - 0.99999972028149386460427194833755*x(1)*x(3); 
dx(3,1) = 1.0000002062427029159152880311012*x(1)*x(2) - 2.666669496024496766040101647377*x(3); 
end