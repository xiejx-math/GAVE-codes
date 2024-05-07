
clear
close all


m=1000;
n=m;
r=n;
dense=0;


if dense
    B=randn(m,n);
else
    B=eye(n);
end

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
if dense
    [xRK,OutRK]=My_RK_GAVE(A,B,b,alpha,opts);
else
    [xRK,OutRK]=My_RK_AVE(A,b,alpha,opts);
end

%% methods of alternating projections
[xMAP,OutMAP]=My_MAP_GAVE(A,B,b,opts);

%% SLA
%epsilon=1;
%[xSLA,OutSLA]=My_SLA_GAVE(A,-B,b,epsilon,opts);

%% generialized Newton method
if dense
    [xGNM,OutGNM]=My_GNM_GAVE(A,B,b,opts);
else
    [xGNM,OutGNM]=My_GNM_AVE(A,b,opts);
end


%% PIM method
if dense
    [xPIM,OutPIM]=My_PIM_GAVE(A,B,b,opts);
else
    [xPIM,OutPIM]=My_PIM_AVE(A,b,opts);
end

%%


fprintf('Iter and CPU: PIM, GNM, MAP, RK, n=%d \n',n)
fprintf(' %4.2f&  %4.2e &  %4.2f &  %4.2e & %4.2f&  %4.2e &  %4.2f &  %4.2e\n',OutPIM.iter,...
    OutPIM.times(end),OutGNM.iter,OutGNM.times(end),OutMAP.iter,OutMAP.times(end),...
    OutRK.iter,OutRK.times(end))




