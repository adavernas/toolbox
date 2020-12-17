function V = matview(A)

V = [A(1,1)-A(2,2) 2*A(1,2)]/(A(1,1)+A(2,2));