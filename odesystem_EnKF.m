function [ f ] = odesystem_EnKF(t, U, G, y, Gamma, I,variance_inflation)
    J = length(U)/(I-1);
    Uhelp = reshape(U,[I-1,J]);
    U = Uhelp;

    Gu = G*U;
    uquer = mean(U,2);
    
    Cuu = (1/J)*U*U'-uquer*uquer';
    f = (Cuu+variance_inflation*1/(t+1)^(1/2)*eye(I-1,I-1))*G'*(Gamma\(y-Gu));
    
    
    fhelp = zeros((I-1)*J,1); %fhelp = f(:);
    for n = 1:J
        fhelp((n-1)*(I-1)+1:n*(I-1),1) = f(:,n).';
    end
    f = fhelp;
end


