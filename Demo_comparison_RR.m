% This Matlab file is used to compare PIM, GNM, MAP, 
% and RABK for the asymmetric ridge regression.
% This file can reproduce Figure 5 in the manuscript

clear
close all

m1=500;
n=200;

mu=500;
lambda=200;

run_time=20;

valuem=m1:m1:10*m1;
sizeP=length(valuem);

RKtimes=zeros(run_time,sizeP);
RKMtimes=zeros(run_time,sizeP);
PIMtimes=zeros(run_time,sizeP);
MAPtimes=zeros(run_time,sizeP);
GNMtimes=zeros(run_time,sizeP);

RKiter=zeros(run_time,sizeP);
RKMiter=zeros(run_time,sizeP);
PIMiter=zeros(run_time,sizeP);
MAPiter=zeros(run_time,sizeP);
GNMiter=zeros(run_time,sizeP);


for jj=1:sizeP
    m=valuem(jj);
    %n=m/2;
    RK_CPU=zeros(run_time,1);
    RKM_CPU=zeros(run_time,1);
    PIM_CPU=zeros(run_time,1);
    MAP_CPU=zeros(run_time,1);
    GNM_CPU=zeros(run_time,1);

    RK_iter=zeros(run_time,1);
    RKM_iter=zeros(run_time,1);
    PIM_iter=zeros(run_time,1);
    MAP_iter=zeros(run_time,1);
    GNM_iter=zeros(run_time,1);
    for ii=1:run_time
        %%
        L=10*rand(m,n)-5*ones(m,n);
        %L=10*randn(m,n);
        c=10*rand(m,1)-5*ones(m,1);
        B=eye(n);
        A=(L'*L+(mu+lambda)*B)/(mu-lambda);
        b=L'*c/(mu-lambda);
        %%
        %opts.xstar=x;
        opts=[];

        %% randomized Kaczmarz method

        alpha=1.0;
        [xRK,OutRK]=My_RK_AVE(A,b,alpha,opts);

        [xRKM,OutRKM]=My_RK_AVE_modified(A,b,alpha,opts);

        %% methods of alternating projections
        [xMAP,OutMAP]=My_MAP_AVE(A,B,b,opts);

        %% generialized Newton method
        [xGNM,OutGNM]=My_GNM_AVE(A,b,opts);

        %% PIM method
        [xPIM,OutPIM]=My_PIM_AVE(A,b,opts);
        %%
        RK_CPU(ii)=OutRK.times(end);
        RKM_CPU(ii)=OutRKM.times(end);
        PIM_CPU(ii)=OutPIM.times(end);
        MAP_CPU(ii)=OutMAP.times(end);
        GNM_CPU(ii)=OutGNM.times(end);

        RK_iter(ii)=OutRK.iter;
        RKM_iter(ii)=OutRKM.iter;
        PIM_iter(ii)=OutPIM.iter;
        MAP_iter(ii)=OutMAP.iter;
        GNM_iter(ii)=OutGNM.iter;

        %%
        RKtimes(ii,jj)=OutRK.times(end);
        RKMtimes(ii,jj)=OutRKM.times(end);
        PIMtimes(ii,jj)=OutPIM.times(end);
        MAPtimes(ii,jj)=OutMAP.times(end);
        GNMtimes(ii,jj)=OutGNM.times(end);

        RKiter(ii,jj)=OutRK.iter/n;
        RKMiter(ii,jj)=OutRKM.iter/n;
        PIMiter(ii,jj)=OutPIM.iter*n;
        MAPiter(ii,jj)=OutMAP.iter*(2*n);
        GNMiter(ii,jj)=OutGNM.iter*n;

    end
    %%
    fprintf(' %d & %4.2f&  %4.2e &  %4.2f &  %4.2e & %4.2f&  %4.2e &  %4.2f &  %4.2e&  %4.2e\n',m,mean(PIM_iter),...
        mean(PIM_CPU),mean(GNM_iter),mean(GNM_CPU),mean(MAP_iter),mean(MAP_CPU),...
        mean(RK_iter),mean(RK_CPU),mean(RKM_CPU))

end


%% plot the CPU time
xlable=valuem;

%%
y1=PIMtimes';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%
y2=GNMtimes';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%
y3=MAPtimes';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);

%%

y4=RKMtimes';
miny4=min(y4');
maxy4=max(y4');
y4q25=quantile(y4,0.25,2);
y4q75=quantile(y4,0.75,2);

%%

figure
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)

h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny4 fliplr(maxy4)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y4q25' fliplr(y4q75')],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)
p1=semilogy( xlable, median(y1'), 'red', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'PIM');
p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 's', 'DisplayName', 'GNM');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
            'LineStyle', '-','Marker', '^', 'DisplayName', 'MAP');
p4=semilogy( xlable, median(y4'), 'magenta', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'x', 'DisplayName', 'RABK');
ylabel('CPU')
xlabel('The values of $m$','Interpreter', 'latex')
%legend('RABK','location', 'best')
legend([p1 p2 p3,p4],{'PIM','GNM','MAP','RABK'},'Interpreter', 'latex','location', 'best')
txt=title(['$n=$ ',num2str(n),', $\bar{\lambda}=$ ',num2str(lambda),', $\bar{\mu}=$',num2str(mu)]);
set(txt, 'Interpreter', 'latex');



%%%%%%%%%%%%%%%%%%%


xlabel_i=valuem;
num_iter_array=xlabel_i';

%% plot errors
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



%%
display_names = {'PIM','GNM','MAP','RABK'};
arrsIter = {PIMiter',GNMiter',MAPiter',RKMiter'};
num_methods = length(arrsIter);
line_colors = {'red','blue','green','magenta'};
minmax_colors = { lightred,lightblue,lightgreen,lightmagenta};
quant_colors = { mediumred,mediumblue,mediumgreen,mediummagenta};

display_legend = true;
max_val_in_plot = 1000000;

%%
[x_arrays_iter, quantiles_iter] =  compute_and_plot4_quantiles_in_logscale(num_iter_array, arrsIter, ...
    num_methods, line_colors, display_names, ...
    minmax_colors, quant_colors, display_legend, max_val_in_plot);
ylabel('Number of full iterations')
xlabel('The values of $m$','Interpreter', 'latex')
txt=title(['$n=$ ',num2str(n),', $\bar{\lambda}=$ ',num2str(lambda),', $\bar{\mu}=$',num2str(mu)]);
set(txt, 'Interpreter', 'latex');


