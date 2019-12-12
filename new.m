clear;


load('reference_data_pde2.mat');



Gamma = eye(K,K)*0.1^2;

T = 10^1;
tspan = [0 T];



w = size(C0,1);
% [V,D]=eig(C0);
% uall = zeros(1,w*w);
% for j = 1:w
%     uall(1,w*(j-1)+1:w*j)=sqrt(D(j,j))*randn*V(:,j);
% end
% save('uall_new2.mat','uall');
load('uall_new.mat');
% 
% uall = 1/100*u0(:);

% a= min(uall);
a = min(utrue)+0.05;
% b = max(uall);
b = max(utrue)-0.05;
% a = min(a1,a2);
% b = max(b1,b2);

% a = -0.3;
% b = 0.1;
T=10^4;

Box = [a,b];
LB = Box(1)*ones(I-1,1);
UB = Box(2)*ones(I-1,1);


funObj = @(w)norm((Gamma^(1/2))\(G*w-y),2)^2;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunctionEvaluation',10^6);
u_KKT = fmincon(funObj,zeros(I-1,1),[],[],[],[],LB,UB,[],options);
% save('u_KKT.mat','u_KKT');
% load('u_KKT.mat');


tic;
J = 15;

    
    u0 = zeros(1,w*J); 
    for j = 1:J
        u0(1,w*(j-1)+1:w*j) = uall(1,w*(j-1)+1:w*j);
    end
    uall = u0;

    variance_inflation = 0.5;
    variance_inflation2 = 0.5;
    alpha = 1/2;
    [t, U, t2, U2, t3, U3] = simulation_new(G,C0,utrue,uall,y,J,tspan,Gamma,I,Box,variance_inflation,variance_inflation2,alpha,u_KKT);
    

time = toc;

save('Estimation.mat','t','U','t2','U2','t3','U3');
save('Information.mat','J','u_KKT','I','K','variance_inflation','variance_inflation2','alpha','Box','T','Gamma','utrue','y','G','C0','time')
