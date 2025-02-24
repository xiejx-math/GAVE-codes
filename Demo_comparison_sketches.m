% This Matlab file is used to compare RIMCS, RIMGS, RIMSRHT, and RABK
% This file can reproduce Figure 1 in the manuscript

clear;
close all;

%%%% size of the matrices
m=2^8;
n=2^7;
kappaA=2;
kappaB=10;
ell=10;

%%%
run_time=20;
Max_length=699;

RABK_CPU=zeros(run_time,Max_length+1);
RABK_error=zeros(run_time,Max_length+1);

RIMUS_CPU=zeros(run_time,Max_length+1);
RIMUS_error=zeros(run_time,Max_length+1);

RIMGS_CPU=zeros(run_time,Max_length+1);
RIMGS_error=zeros(run_time,Max_length+1);

RIMCS_CPU=zeros(run_time,Max_length+1);
RIMCS_error=zeros(run_time,Max_length+1);

RIMSRHT_CPU=zeros(run_time,Max_length+1);
RIMSRHT_error=zeros(run_time,Max_length+1);
%%%% run and store the numerical results

for ii=1:run_time
    [A,B]=getAB(m,n,2,kappaA,1,kappaB);
    x=randn(n,1);
    b=A*x-B*abs(x);

    %%%%
    opts.xstar=x;
    opts.Max_iter=Max_length;
    opts.TOL=eps^2;

    %%% randomized averaged block Kaczmarz method

    alpha=1.0;
    [xRABK,OutRABK]=My_RABK1_GAVE(A,B,b,alpha,ell,opts);

    %%% randomized iterative method with Guassian sketching
    [xRIMGS,OutRIMGS]=My_RIMGS_GAVE(A,B,b,alpha,ell,opts);

    %%% randomized iterative method with Countsketch
    [xRIMCS,OutRIMCS]=My_RIMCS_GAVE(A,B,b,alpha,ell,opts);

    %%% randomized iterative method with SRHT
    [xRIMSRHT,OutRIMSRHT]=My_RIMSRHT_GAVE(A,B,b,alpha,ell,opts);

    %
    RABK_error(ii,1:length(OutRABK.error))=OutRABK.error;
    RABK_CPU(ii,1:length(OutRABK.times))=OutRABK.times;

    RIMGS_error(ii,1:length(OutRIMGS.error))=OutRIMGS.error;
    RIMGS_CPU(ii,1:length(OutRIMGS.times))=OutRIMGS.times;

    RIMCS_error(ii,1:length(OutRIMCS.error))=OutRIMCS.error;
    RIMCS_CPU(ii,1:length(OutRIMCS.times))=OutRIMCS.times;

    RIMSRHT_error(ii,1:length(OutRIMSRHT.error))=OutRIMSRHT.error;
    RIMSRHT_CPU(ii,1:length(OutRIMSRHT.times))=OutRIMSRHT.times;

    fprintf('Done, ii=%d\n',ii)
end


%%%
xlable=1:(Max_length+1);

%%%
y1=RABK_error';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%%
y2=RIMCS_error';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);


%%%
y3=RIMGS_error';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);



%%%
y4=RIMSRHT_error';
miny4=min(y4');
maxy4=max(y4');
y4q25=quantile(y4,0.25,2);
y4q75=quantile(y4,0.75,2);

%%%
figure
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)

h = fill([xlable  fliplr(xlable)], [miny4 fliplr(maxy4)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y4q25' fliplr(y4q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)

h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)
%%%
p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
    'LineStyle', '--', 'DisplayName', 'RIMCS');
p3=semilogy( xlable, median(y3'), 'red', 'LineWidth', 1,...
    'LineStyle', '-.', 'DisplayName', 'RIMGS');
p4=semilogy( xlable, median(y4'), 'green', 'LineWidth', 1,...
    'LineStyle', '-', 'DisplayName', 'RIMSRHT');

p1=semilogy( xlable, median(y1'), 'magenta', 'LineWidth', 1,...
    'LineStyle', ':', 'DisplayName', 'RABK');
set(gca, 'YScale', 'log')
ylim([10^(-12), 1])
xlim([0, 650])
ylabel('RSE','Interpreter', 'latex')
xlabel('Number of iterations','Interpreter', 'latex')
legend([ p2 p3 p4 p1],{'RIMCS','RIMGS','RIMSRHT','RABK'},'Interpreter', 'latex','location', 'best')
txt=title(['$m=$ ',num2str(m),', $n=$ ',num2str(n),', $\kappa_A=$ ',num2str(kappaA),', $\kappa_B=$ ',num2str(kappaB)]);
set(txt, 'Interpreter', 'latex');





%%%

figure
h = fill([median(RIMCS_CPU)  fliplr(median(RIMCS_CPU))], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([median(RIMCS_CPU)  fliplr(median(RIMCS_CPU))], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([median(RIMGS_CPU)  fliplr(median(RIMGS_CPU) )], [miny3 fliplr(maxy3)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([median(RIMGS_CPU)  fliplr(median(RIMGS_CPU) )], [y3q25' fliplr(y3q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)

h = fill([median(RIMSRHT_CPU)  fliplr(median(RIMSRHT_CPU))], [miny4 fliplr(maxy4)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([median(RIMSRHT_CPU)  fliplr(median(RIMSRHT_CPU))], [y4q25' fliplr(y4q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)

h = fill([median(RABK_CPU)  fliplr(median(RABK_CPU) )], [miny1 fliplr(maxy1)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([median(RABK_CPU)  fliplr(median(RABK_CPU) )], [y1q25' fliplr(y1q75')],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)
%%%
p2=semilogy( median(RIMCS_CPU), median(y2'), 'blue', 'LineWidth', 1,...
    'LineStyle', '--', 'DisplayName', 'RIMCS');
p3=semilogy( median(RIMGS_CPU), median(y3'), 'red', 'LineWidth', 1,...
    'LineStyle', '-.', 'DisplayName', 'RIMGS');
p4=semilogy( median(RIMSRHT_CPU), median(y4'), 'green', 'LineWidth', 1,...
    'LineStyle', '-', 'DisplayName', 'RIMSRHT');

p1=semilogy( median(RABK_CPU), median(y1'), 'magenta', 'LineWidth', 1,...
    'LineStyle', ':', 'DisplayName', 'RABK');
set(gca, 'YScale', 'log')
xlim([0,2.2]);
ylim([10^(-12), 1]);
ylabel('RSE','Interpreter', 'latex')
xlabel('CPU','Interpreter', 'latex')
legend([ p2 p3 p4 p1],{'RIMCS','RIMGS','RIMSRHT','RABK'},'Interpreter', 'latex','location', 'best')
txt=title(['$m=$ ',num2str(m),', $n=$ ',num2str(n),', $\kappa_A=$ ',num2str(kappaA),', $\kappa_B=$ ',num2str(kappaB)]);
set(txt, 'Interpreter', 'latex');
axes('position',[0.6,0.45,0.28,0.28]);
h = fill([median(RIMCS_CPU)  fliplr(median(RIMCS_CPU))], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([median(RIMCS_CPU)  fliplr(median(RIMCS_CPU))], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([median(RIMGS_CPU)  fliplr(median(RIMGS_CPU) )], [miny3 fliplr(maxy3)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([median(RIMGS_CPU)  fliplr(median(RIMGS_CPU) )], [y3q25' fliplr(y3q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)

h = fill([median(RABK_CPU)  fliplr(median(RABK_CPU) )], [miny1 fliplr(maxy1)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([median(RABK_CPU)  fliplr(median(RABK_CPU) )], [y1q25' fliplr(y1q75')],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)
set(h,'facealpha', .1)
p2=semilogy( median(RIMCS_CPU), median(y2'), 'blue', 'LineWidth', 1,...
    'LineStyle', '--', 'DisplayName', 'RIMCS');
p3=semilogy( median(RIMGS_CPU), median(y3'), 'red', 'LineWidth', 1,...
    'LineStyle', '-.', 'DisplayName', 'RIMGS');
p1=semilogy( median(RABK_CPU), median(y1'), 'magenta', 'LineWidth', 1,...
    'LineStyle', ':', 'DisplayName', 'RABK');
xlim([0,0.28]);
ylim([10^(-12), 1]);
set(gca, 'YScale', 'log')

