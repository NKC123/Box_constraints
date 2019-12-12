%%% Evaluation
clear all;

% load data
load('Estimation.mat');
load('Information.mat');

w = size(C0,1);
xx = linspace(0,pi,w);

% set figure(1) for estimation in parameter space
fig1 = figure(1);    
clf(fig1)
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25]);

% plot the true parameter
plot(xx,utrue,'g','DisplayName','truth');hold on

% compute the projected true parameter
Putrue = utrue;
Putrue(Putrue<Box(1))= Box(1);
Putrue(Putrue>Box(2))= Box(2);

% plot the projected true parameter
plot(xx,Putrue,'black','DisplayName','projected truth');hold on
plot(xx,u_KKT,'magenta','DisplayName','$u^*$','LineWidth',2);hold on
xlabel('x')


% set figure(2) for estimation in observation space
% plot the true observations
fig2 = figure(2);
clf(fig2)
set(fig2, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25]);
plot(linspace(0,pi,length([0;y;0])),[0;y;0],'go','DisplayName','observations','LineWidth',2);hold on
xlabel('x')


% compute the estimates of the run simulations
Uhelp = reshape(U(end,:),[I-1,J]);
Uest = mean(Uhelp,2);
Uhelp2 = reshape(U2(end,:),[I-1,J]);
Uest2 = mean(Uhelp2,2);
Uhelp3 = reshape(U3(end,:),[I-1,J]);
Uest3 = mean(Uhelp3,2);


% plot the EnKF estimates in the parameter space
figure(1)
plot(xx,Uest3,'-','DisplayName',['transformed EnKF, J=' num2str(J)],'LineWidth',2);hold on
plot(xx,Uest,'-.','DisplayName',['EnKF without BC, J=' num2str(J)],'LineWidth',2);hold on
plot(xx,Uest2,'--','DisplayName',['projected EnKF, J=' num2str(J)],'LineWidth',2);hold on
l = legend('show','Location','northeast');
set(l,'Interpreter','latex','FontSize',15);

% plot the EnKF estimates in the observation space
figure(2)
plot(linspace(0,pi,length([0;y;0])),[0;G*Uest3;0],'+','DisplayName',['transformed EnKF, J=' num2str(J)],'LineWidth',2);hold on
plot(linspace(0,pi,length([0;y;0])),[0;G*Uest;0],'.','DisplayName',['EnKF without BC, J=' num2str(J)],'LineWidth',2);hold on
plot(linspace(0,pi,length([0;y;0])),[0;G*Uest2;0],'.','DisplayName',['projected EnKF, J=' num2str(J)],'LineWidth',2);hold on
l = legend('show','Location','northeast');
set(l,'Interpreter','latex','FontSize',15);



%%%% compute the different error terms e and r:


%%% EnKF without BC
% ensemble spread in the parameter space
e = zeros(1,length(t));
% residuals in the parameter space
r = zeros(1,length(t));
% ensemble spread in the observation space
Ae = zeros(1,length(t));
% residuals in the observation space
Ar = zeros(1,length(t));
% residuals of the cost function
Phi_diff =zeros(1,length(t));

%%% projected EnKF
% ensemble spread in the parameter space
e2 = zeros(1,length(t2));
% residuals in the parameter space
r2 = zeros(1,length(t2));
% ensemble spread in the observation space
Ae2 = zeros(1,length(t2));
% residuals in the observation space
Ar2 = zeros(1,length(t2));
% difference to the KKT-point (not necessarily unique) in the parameter
% space
Pr2 = zeros(1,length(t2));
% difference to the KKT-point (not necessarily unique) in the observation
% space
PAr2 = zeros(1,length(t2));
% residuals of the cost function
Phi_diff2 =zeros(1,length(t2));

%%% Transformed EnKF
% ensemble spread in the parameter space
e3 = zeros(1,length(t3));
% residuals in the parameter space
r3 = zeros(1,length(t3));
% ensemble spread in the observation space
Ae3 = zeros(1,length(t3));
% residuals in the observation spac
Ar3 = zeros(1,length(t3));
% difference to the KKT-point (not necessarily unique) in the parameter
% space
Pr3 = zeros(1,length(t3));
% difference to the KKT-point (not necessarily unique) in the observation
% space
PAr3 = zeros(1,length(t3));
% residuals of the cost function
Phi_diff3 =zeros(1,length(t3));

% cost function value of the KKT-point
r_KKT = u_KKT-utrue;
Phi_KKT = norm((Gamma^(1/2))\(G*r_KKT),2)^2;

%%% does the mean converge?
% Cauchy-differences
conv_mean1=zeros(1,length(t)-1);
conv_mean2=zeros(1,length(t2)-1);
conv_mean3=zeros(1,length(t3)-1);

% Computation of the EnKF without BC
for k = 1:length(t)
    Uhelp = reshape(U(k,:),[I-1,J]);
    Uquer = mean(Uhelp,2);
    edt = Uhelp-Uquer;
    rdt = Uhelp-utrue;
    Aedt = (Gamma^(1/2))\(G*edt);
    Ardt = (Gamma^(1/2))\(G*rdt);
    e(1,k) = norm(edt(:),2)^2;
    r(1,k) = norm(rdt(:),2)^2;
    Ae(1,k) = norm(Aedt(:),2)^2;
    Ar(1,k) = norm(Ardt(:),2)^2;
    Phi_u = norm(Gamma^(1/2)\(G*rdt),2)^2;
    Phi_diff(1,k) = (1/J*Phi_u-Phi_KKT)^2;
    if k>=2
        Uhelp2 = reshape(U(k-1,:),[I-1,J]);
        Uquer2 = mean(Uhelp2,2);
        conv_mean1(1,k-1)=norm(Uquer-Uquer2,2)^2;
    end
end

% Computation of the projected EnKF
for k=1:length(t2) 
    Uhelp2 = reshape(U2(k,:),[I-1,J]);
    Uquer2 = mean(Uhelp2,2);
    edt2 = Uhelp2-Uquer2;
    rdt2 = Uhelp2-utrue;
    Prdt2 = Uhelp2-u_KKT;
    Aedt2 = (Gamma^(1/2))\(G*edt2);
    Ardt2 = (Gamma^(1/2))\(G*rdt2);
    PArdt2 = (Gamma^(1/2))\(G*Prdt2);
    e2(1,k) = norm(edt2(:),2)^2;
    r2(1,k) = norm(rdt2(:),2)^2;
    Pr2(1,k) = norm(Prdt2(:),2)^2;
    Ae2(1,k) = norm(Aedt2(:),2)^2;
    Ar2(1,k) = norm(Ardt2(:),2)^2;
    PAr2(1,k) = norm(PArdt2(:),2)^2;
    Phi_u = norm(Gamma^(1/2)\(G*rdt2),2)^2;
    Phi_diff2(1,k) = (1/J*Phi_u-Phi_KKT)^2;
    if k>=2
        Uhelp2 = reshape(U2(k-1,:),[I-1,J]);
        Uquer2 = mean(Uhelp2,2);
        conv_mean2(1,k-1)=norm(Uquer-Uquer2,2)^2;
    end
end

for k=1:length(t3) 
    Uhelp3 = reshape(U3(k,:),[I-1,J]);
    Uquer3 = mean(Uhelp3,2);
    edt3 = Uhelp3-Uquer3;
    rdt3 = Uhelp3-utrue;
    Prdt3 = Uhelp3-u_KKT;
    Aedt3 = (Gamma^(1/2))\(G*edt3);
    Ardt3 = (Gamma^(1/2))\(G*rdt3);
    PArdt3 = (Gamma^(1/2))\(G*Prdt3);
    e3(1,k) = norm(edt3(:),2)^2;
    r3(1,k) = norm(rdt3(:),2)^2;
    Pr3(1,k) = norm(Prdt3(:),2)^2;
    Ae3(1,k) = norm(Aedt3(:),2)^2;
    Ar3(1,k) = norm(Ardt3(:),2)^2;
    PAr3(1,k) = norm(PArdt3(:),2)^2;
    Phi_u = norm((Gamma^(1/2))\(G*rdt3),2)^2;
    Phi_diff3(1,k) = (1/J*Phi_u-Phi_KKT)^2;
    if k>=2
        Uhelp2 = reshape(U3(k-1,:),[I-1,J]);
        Uquer2 = mean(Uhelp2,2);
        conv_mean3(1,k-1)=norm(Uquer-Uquer2,2)^2;
    end
end

% set figure(3) for the ensemble spread in parameter space
fig3 = figure(3);
clf(fig3)
set(fig3, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25]);
formatlegend = '$\\frac{1}{J}\\sum |e_t^{(j)}|^2$ transformed EnKF';
str = sprintf(formatlegend);
loglog(t3,e3,'-','DisplayName',str,'LineWidth',2); hold on
formatlegend = '$\\frac{1}{J}\\sum |e_t^{(j)}|^2$ EnKF without BC';
str = sprintf(formatlegend);
loglog(t,e,'--','DisplayName',str,'LineWidth',2); hold on
formatlegend = '$\\frac{1}{J}\\sum |e_t^{(j)}|^2$ projected EnKF';
str = sprintf(formatlegend);
loglog(t2,e2,'-.','DisplayName',str,'LineWidth',2); hold on
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',15);
title('Parameter space')

% set figure(4) for the residuals in parameter space
fig4 = figure(4);
clf(fig4)
set(fig4, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25]);
formatlegend = '$\\frac{1}{J}\\sum |r_t^{(j)}|^2$ transformed EnKF';
str = sprintf(formatlegend);
loglog(t3,r3,'-','DisplayName',str,'LineWidth',2);hold on
formatlegend = '$\\frac{1}{J}\\sum |r_t^{(j)}|^2$ EnKF without BC';
str = sprintf(formatlegend);
loglog(t,r,'--','DisplayName',str,'LineWidth',2);hold on
formatlegend = '$\\frac{1}{J}\\sum |r_t^{(j)}|^2$ projected EnKF';
str = sprintf(formatlegend);
loglog(t2,r2,'-.','DisplayName',str,'LineWidth',2);hold on
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',15);
title('Parameter space')

% set figure(5) for the ensemble spread in observation space
fig5 = figure(5);
clf(fig5)
set(fig5, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25]);
formatlegend = '$\\frac{1}{J}\\sum |Ae_t^{(j)}|_{\\Gamma}^2$ transformed EnKF';
str = sprintf(formatlegend);
loglog(t3,Ae3,'-','DisplayName',str,'LineWidth',2);hold on
formatlegend = '$\\frac{1}{J}\\sum |Ae_t^{(j)}|_{\\Gamma}^2$ EnKF without BC';
str = sprintf(formatlegend);
loglog(t,Ae,'--','DisplayName',str,'LineWidth',2);hold on
formatlegend = '$\\frac{1}{J}\\sum |Ae_t^{(j)}|_{\\Gamma}^2$ projected EnKF';
str = sprintf(formatlegend);
loglog(t2,Ae2,'-.','DisplayName',str,'LineWidth',2);hold on
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',15);
title('Observation space')

% set figure(6) for the residuals in observation space
fig6 = figure(6);
clf(fig6)
set(fig6, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25]);
formatlegend = '$\\frac{1}{J}\\sum |Ar_t^{(j)}|_{\\Gamma}^2$ transformed EnKF';
str = sprintf(formatlegend);
loglog(t3,Ar3,'-','DisplayName',str,'LineWidth',2);hold on
formatlegend = '$\\frac{1}{J}\\sum |Ar_t^{(j)}|_{\\Gamma}^2$ EnKF without BC';
str = sprintf(formatlegend);
loglog(t,Ar,'--','DisplayName',str,'LineWidth',2);hold on
formatlegend = '$\\frac{1}{J}\\sum |Ar_t^{(j)}|_{\\Gamma}^2$ projected EnKF';
str = sprintf(formatlegend);
loglog(t2,Ar2,'-.','DisplayName',str,'LineWidth',2);hold on
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',15);
title('Observation space')

% set figure(7) for the difference to the KKT-point in the parameter space
fig7 = figure(7);
clf(7)
set(fig7, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25]);
formatlegend = '$\\frac{1}{J}\\sum |u_t^{(j)}-u^*|^2$ transformed EnKF';
str = sprintf(formatlegend);
loglog(t3,Pr3,'-','DisplayName',str,'LineWidth',2);hold on
formatlegend = '$\\frac{1}{J}\\sum |u_t^{(j)}-u^*|^2$ projected EnKF';
str = sprintf(formatlegend);
loglog(t2,Pr2,'-.','DisplayName',str,'LineWidth',2);hold on
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',15);
title('Parameter space')

% set figure(8) for the difference to the KKT-point in the observation space
fig8 = figure(8);
clf(8)
set(fig8, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25]);
formatlegend = '$\\frac{1}{J}\\sum |A(u_t^{(j)}-u^*)|_{\\Gamma}^2$ transformed EnKF';
str = sprintf(formatlegend);
loglog(t3,PAr3,'-','DisplayName',str,'LineWidth',2);hold on
formatlegend = '$\\frac{1}{J}\\sum |A(u_t^{(j)}-u^*)|_{\\Gamma}^2$ projected EnKF';
str = sprintf(formatlegend);
loglog(t2,PAr2,'-.','DisplayName',str,'LineWidth',2);hold on
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',15);
title('Observation space')

% set figure(9) for the difference of the costfunction to the KKT-point
fig9 = figure(9);
clf(9)
set(fig9, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25]);
formatlegend = '$\\frac{1}{J}\\sum (\\Phi(u_t^{(j)})-\\Phi(u^*))^2$ transformed EnKF';
str = sprintf(formatlegend);
loglog(t3,Phi_diff3,'-','DisplayName',str,'LineWidth',2);hold on
formatlegend = '$\\frac{1}{J}\\sum (\\Phi(u_t^{(j)})-\\Phi(u^*))^2$ EnKF without BC';
str = sprintf(formatlegend);
loglog(t,Phi_diff,'-.','DisplayName',str,'LineWidth',2);hold on
formatlegend = '$\\frac{1}{J}\\sum (\\Phi(u_t^{(j)})-\\Phi(u^*))^2$ projected EnKF';
str = sprintf(formatlegend);
loglog(t2,Phi_diff2,'-.','DisplayName',str,'LineWidth',2);hold on
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',15);

% set figure(10) to test convergence of the mean
fig10 = figure(10);
clf(10)
set(fig10, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25]);
formatlegend = 'Convergence $\\bar{u}_t$ transformed EnKF';
str = sprintf(formatlegend);
loglog(t3(2:end),conv_mean3,'-','DisplayName',str,'LineWidth',2);hold on
formatlegend = 'Convergence $\\bar{u}_t$ EnKF without BC';
str = sprintf(formatlegend);
loglog(t(2:end),conv_mean1,':','DisplayName',str,'LineWidth',2);hold on
formatlegend = 'Convergence $\\bar{u}_t$ projected EnKF';
str = sprintf(formatlegend);
loglog(t2(2:end),conv_mean2,'-.','DisplayName',str,'LineWidth',2);hold on
l = legend('show','Location','southwest');
set(l,'Interpreter','latex','FontSize',15);
