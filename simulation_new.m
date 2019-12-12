function [t,U,t2,U2,t3,U3 ] = simulation_new(G,C0,utrue,u0,y,J,tspan,Gamma,I,Box,variance_inflation, variance_inflation2,alpha,u_KKT)
    
    w = size(C0,1);
    xx = linspace(0,pi,w);
    
    % time
    u2 = u0;
    u2(u2<Box(1)) = Box(1);
    u2(u2>Box(2)) = Box(2);

    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [t, U] = ode45(@(t,U) odesystem_EnKF(t,U,G,y,Gamma,I,variance_inflation2), tspan, u0);
    [t2, U2] = ode45(@(t2,U2) odesystem_EnKF_BC_correction(t2,U2,G,y,Gamma,I,Box,variance_inflation2), tspan, u2);
    [t3, U3] = ode45(@(t3,U3) odesystem_new(t3,U3,G,y,Gamma,I,Box,variance_inflation,alpha), tspan, u2,options);
end


