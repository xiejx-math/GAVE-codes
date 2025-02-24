function [A,B]=getAB(m,n,sigma_A_min,kappa_A,sigma_B_max,kappa_B)

%%%
sigma_A = linspace(sigma_A_min, kappa_A*sigma_A_min, n);
A = diag(sigma_A); 
%%%
sigma_B = linspace(sigma_B_max/kappa_B, sigma_B_max, n); 
B = diag(sigma_B); 
%%%
[U_A,~]=qr(randn(m,n),0);
[V_A,~]=qr(randn(n,n),0);
[U_B,~]=qr(randn(m,n),0);
[V_B,~]=qr(randn(n,n),0);
%%%
A = U_A * A * V_A'; 
B = U_B * B * V_B'; 

end