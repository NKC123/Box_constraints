%%% Karhunen Loeve expansion (Gaussian random field) %%%
%%% Based on the Matern covariance family %%%
%%% Input the parameters alpha and tau, NOTE: alpha must be chosen > d/2 %%%
%%% To view the random field after L is computed, input:
%%% surf(idst2(reshape(L,N,N))),'edgecolor','none');%%%
%%% KLE is computed through the eigenvalues and eigenvectors of the covariance operator %%%


function L = kle(alpha,tau,N)
	
	% Random variables in KL expansion
    xi = normrnd(0,1,N);
   % xi = speye(N);
	%xi=1;
	% Define the (square root of) eigenvalues of the covariance operator
	[K1,K2] = meshgrid(0:N-1,0:N-1);
	coef = (pi^2*(K1.^2+K2.^2) + tau^2).^(-alpha/2);	
	
	% Construct the KL coefficients
	L = N*coef.*xi;
    L(1,1) = 0;
	
    L = reshape(L,N^2,1);
    
end


