function [ f ] = odesystem_new( t,U, G, y, Gamma, I,Box,variance_inflation,alpha )
    
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
    Cuu = (1/J)*U*U'-uquer*uquer';
%     Cup = Cuu*G';
%     f = zeros(I-1,J);
%     for n=1:J
%         Cuu_n = Cuu+variance_inflation*eye(I-1,I-1);
%         K1 = (Uhelp(:,n)<=Box(1));
%         K2 = (Uhelp(:,n)>=Box(2));
%         for l1 = K1
%             del = 1:(I-1);
%             del(l1) = [];
%             Cuu_n(del,l1) = 0;
%             Cuu_n(l1,del) = 0;
%         end
%         for l1 = K2
%             del = 1:(I-1);
%             del(l1) = [];
%             Cuu_n(del,l1) = 0;
%             Cuu_n(l1,del) = 0;
%         end
%         f(:,n) = Cuu_n*G'*(Gamma\(y-G*U(:,n)));
%     end
    Cuu_n = Cuu+variance_inflation*eye(I-1,I-1);%1/(t+1)^(alpha)*variance_inflation*eye(I-1,I-1);
    for n=1:J
        
%         K1 = (Uhelp(:,n)<=Box(1));
%         K2 = (Uhelp(:,n)>=Box(2));
%         Test = G'*(Gamma\(G*Uhelp(:,n)-y));
%         K3 = (Test(K1)<0);
%         K4 = (Test(K2)>0);
%         K1(K3) = [];
%         K2(K4) = [];

        K1 = (Uhelp(:,n)<=Box(1));
        K2 = (Uhelp(:,n)>=Box(2));
        Test = G'*(Gamma\(G*Uhelp(:,n)-y));
        K3 = (Test<=0);
        K4 = (Test>=0);
        K1(K1==K3) = 0;
        K2(K2==K4) = 0;
          
        for l1 = K1
            del = 1:(I-1);
            del(l1) = [];
            Cuu_n(del,l1) = 0;
            Cuu_n(l1,del) = 0;
        end
        for l1 = K2
            del = 1:(I-1);
            del(l1) = [];
            Cuu_n(del,l1) = 0;
            Cuu_n(l1,del) = 0;
        end
    end
    f = Cuu_n*G'*(Gamma\(y-Gu));
    
%     f = Cup*(Gamma\(y-1/2*Gu-1/2*Gquer));
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