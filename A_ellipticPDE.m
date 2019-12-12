function [ A ] = A_ellipticPDE( N )
    h = pi/(N+1);
    A = zeros(N,N);
    A(1,1:2) = [2+h^2, -1];
    A(N,(N-1):N) = [-1, 2+h^2];
    for k = 2:(N-1)
        A(k,(k-1):(k+1)) = [-1, 2+h^2, -1];
    end
    A = h^(-2)*A;
end

