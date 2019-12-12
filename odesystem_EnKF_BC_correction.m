function [ f ] = odesystem_EnKF_BC_correction( t,U, G, y, Gamma, I,Box,variance_inflation )
    
    I1 = (U<=Box(1));
    I2 = (U>=Box(2));

    U(I1) = Box(1);
    U(I2) = Box(2);
    
% %     b1 = Box(1,:);
% %     b2 = Box(2,:);
% %     Ind1 = (U<b1);
% %     U(Ind1) = b1(Ind1);
% %     Ind2 = (U>b2);
% %     U(Ind2) = b2(Ind2);

    J = length(U)/(I-1);
    Uhelp = reshape(U,[I-1,J]);
    U = Uhelp;

    Gu = G*U;
%     Gquer = mean(Gu,2);
    uquer = mean(U,2);
    
%     Cup = (1/J)*U*Gu'-uquer*Gquer';
%     f = Cup*(Gamma\(y-Gu));
    Cuu = (1/J)*U*U'-uquer*uquer';
    f = (Cuu+variance_inflation*eye(I-1,I-1))*G'*(Gamma\(y-Gu));%*1/(t+1)^(1/2)*eye(I-1,I-1))*G'*(Gamma\(y-Gu));
    
    for i = I1
        f(i) = max(0,f(i));
    end
    for j = I2
        f(j) = min(0,f(j));
    end
    
    fhelp = zeros((I-1)*J,1); %fhelp = f(:);
    for n = 1:J
        fhelp((n-1)*(I-1)+1:n*(I-1),1) = f(:,n).';
    end
    f = fhelp;
end

