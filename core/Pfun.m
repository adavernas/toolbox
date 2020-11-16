function k = Pfun(N,P,H)

% n is normal
% p is point on plane
% h is point to project
% k is projection

a = N(1);
b = N(2);
c = N(3);

d = P(1);
e = P(2);
f = P(3);

x = H(1);
y = H(2);
z = H(3);

t = (a*d - a*x + b*e - b*y + c*f - c*z)/(a^2+b^2+c^2);

k = [x+t*a; y+t*b; z+t*c];