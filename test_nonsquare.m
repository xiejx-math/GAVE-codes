
clear;
close all;

%% size of the matrices
m=2000;
n=1000;
r=n;


%% run and store the numerical results

B=randn(m,n);
[U,~]=qr(randn(m,r),0);
[V,~]=qr(randn(n,r),0);
sigmaA=norm(B)*eye(r)+diag(rand(r,1));
A=U*sigmaA*V';
clear U V sigmaA
x=randn(n,1);
b=A*x-B*abs(x);

%%
opts.xstar=x;

%% randomized Kaczmarz method

alpha=1.0;
[xRK,OutRK]=My_RK_GAVE(A,B,b,alpha,opts);


%% methods of alternating projections
[xMAP,OutMAP]=My_MAP_GAVE(A,B,b,opts);

%% SLA
epsilon=1;
[xSLA,OutSLA]=My_SLA_GAVE(A,-B,b,epsilon,opts);
%% store the computational results


%% print the result at each step
fprintf('Iter and CPU: SLA, MAP, RK; m=%d, n=%d \n',m,n)
fprintf(' %4.2f &  %4.2f & %4.2f&  %4.3f &  %4.2f &  %4.3f\n',OutSLA.iter,...
    OutSLA.times(end),OutMAP.iter,OutMAP.times(end),...
    OutRK.iter,OutRK.times(end))





