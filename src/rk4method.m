function [ W ] = rk4method( x_init, f, h, T)
% This code was created by Joel Rosenfeld in 2021 to accompany his YouTube
% channel ThatMaththing (http://www.thatmaththing.com/
% If you use this code for a project, please credit Joel A. Rosenfeld, and
% link his YouTube channel and professional website,
% http://www.thelearningdock.org/
%RK4 Method


% f should be a function that returns column vectors
% x_init should be a column vector


W = [x_init];

for i=1:T/h-1
    try
        V1 = f(W(:,i));
        V2 = f(W(:,i) + 1/2*h*V1);
        V3 = f(W(:,i) + 1/2*h*V2);
        V4 = f(W(:,i) + h*V3);
    catch
        V1 = f(0,W(:,i));
        V2 = f(0,W(:,i) + 1/2*h*V1);
        V3 = f(0,W(:,i) + 1/2*h*V2);
        V4 = f(0,W(:,i) + h*V3);
    end
    W = [W,W(:,i) + h/6*(V1+2*V2+2*V3+V4)];
end


end

