% This Matlab file is used to compare SLA, MAP, and RABK
% This file can reproduce Figure 4 in the manuscript

clear;
close all;

%%% size of the matrices
m1=200;
n=200;
kappaA=2;
kappaB=2;
%%% set some paramters
ss=10;
run_time=20;
sizeP=ss;
RKtimes=zeros(run_time,sizeP);
MAPtimes=zeros(run_time,sizeP);
SLAtimes=zeros(run_time,sizeP);

%%% run and store the numerical results
for ii=1:ss
    m=m1+m1*(ii-1);
    RK_CPU=zeros(1,run_time);
    MAP_CPU=zeros(1,run_time);
    SLA_CPU=zeros(1,run_time);
    RK_Iter=zeros(1,run_time);
    MAP_Iter=zeros(1,run_time);
    SLA_Iter=zeros(1,run_time);
    for jj=1:run_time
        [A,B]=getAB(m,n,2,kappaA,1,kappaB);
        x=randn(n,1);
        b=A*x-B*abs(x);
        opts.xstar=x;

        %%% randomized  average block Kaczmarz with p=1 method

        alpha=1.0;
        [xRK,OutRK]=My_RK_GAVE(A,B,b,alpha,opts);

        %%% methods of alternating projections
        [xMAP,OutMAP]=My_MAP_GAVE(A,B,b,opts);

        %%% SLA
        epsilon=1;
        [xSLA,OutSLA]=My_SLA_GAVE(A,-B,b,epsilon,opts);

        %%% store the computational results
        RK_CPU(jj)=OutRK.times(end);
        MAP_CPU(jj)=OutMAP.times(end);
        SLA_CPU(jj)=OutSLA.times(end);

        RK_Iter(jj)=OutRK.iter;
        MAP_Iter(jj)=OutMAP.iter;
        SLA_Iter(jj)=OutSLA.iter;
       
        RKtimes(jj,ii)=OutRK.times(end);
        MAPtimes(jj,ii)=OutMAP.times(end);
        SLAtimes(jj,ii)=OutSLA.times(end);
    end
    %%% print the result at each step
    fprintf('Iter and CPU: SLA, MAP, RK; m=%d, n=%d \n',m,n)
    fprintf(' %4.2f &  %4.2f & %4.2f&  %4.3f &  %4.2f &  %4.3f\n',mean(SLA_Iter),...
        mean(SLA_CPU),mean(MAP_Iter),mean(MAP_CPU),...
        mean(RK_Iter),mean(RK_CPU))
end





%%% plot the CPU time

xlabel_i=m1*[1:ss];
num_iter_array=xlabel_i';

%%% set parameters
lightgray =   [0.8 0.8 0.8];
mediumgray =  [0.6 0.6 0.6];
lightred =    [1 0.9 0.9];
mediumred =   [1 0.6 0.6];
lightgreen =  [0.9 1 0.9];
mediumgreen = [0.6 1 0.6];
lightblue =   [0.9 0.9 1];
mediumblue =  [0.6 0.6 1];
lightmagenta =   [1 0.9 1];
mediummagenta =  [1 0.6 1];

%%%
display_names = {'SLA','MAP','RABK'};
arrsIter = {SLAtimes',MAPtimes',RKtimes'};
num_methods = length(arrsIter);
line_colors = {'black','green','magenta'};
minmax_colors = { lightgray, lightgreen,lightmagenta};
quant_colors = { mediumgray,mediumgreen,mediummagenta};

display_legend = true;
max_val_in_plot = 1000;

%%%
[x_arrays_iter, quantiles_iter] =  compute_and_plot_RR_quantiles_in_logscale(num_iter_array, arrsIter, ...
    num_methods, line_colors, display_names, ...
    minmax_colors, quant_colors, display_legend, max_val_in_plot);
ylabel('CPU')
xlabel('Number of rows $(m)$','Interpreter', 'latex')
txt=title(['$n=$', num2str(n),', $\kappa_A=$',num2str(kappaA),', $\kappa_B=$ ',num2str(kappaB)]);
set(txt, 'Interpreter', 'latex');








