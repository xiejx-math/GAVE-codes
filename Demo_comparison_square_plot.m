% This Matlab file is used to compare PIM, GNM, MAP,
% and RABK with square matrices.
% This file can reproduce Figure 3 in the manuscript

clear
close all

m1=1000;
dense=1;
run_time=20;
Numbn=10;

kappaA=2;
kappaB=1;

CPUdata=zeros(4,Numbn);
RK_CPU=zeros(run_time,Numbn);
PIM_CPU=zeros(run_time,Numbn);
MAP_CPU=zeros(run_time,Numbn);
GNM_CPU=zeros(run_time,Numbn);

RK_iter=zeros(run_time,Numbn);
PIM_iter=zeros(run_time,Numbn);
MAP_iter=zeros(run_time,Numbn);
GNM_iter=zeros(run_time,Numbn);

for jj=1:Numbn
    m=m1*jj;
    n=m;
    for ii=1:run_time
        [A,B]=getAB(m,n,2,kappaA,1,kappaB);
        x=randn(n,1);
        b=A*x-B*abs(x);
        opts.xstar=x;

        %%% randomized average block Kaczmarz with p=1 method
        alpha=1.0;
        [xRK,OutRK]=My_RK_GAVE(A,B,b,alpha,opts);

        %%% methods of alternating projections
        [xMAP,OutMAP]=My_MAP_GAVE(A,B,b,opts);

        %%% generialized Newton method
        [xGNM,OutGNM]=My_GNM_GAVE(A,B,b,opts);

        %%% PIM method
        if dense
            [xPIM,OutPIM]=My_PIM_GAVE(A,B,b,opts);
        else
            [xPIM,OutPIM]=My_PIM_AVE(A,b,opts);
        end

        %%%
        RK_CPU(ii,jj)=OutRK.times(end);
        PIM_CPU(ii,jj)=OutPIM.times(end);
        MAP_CPU(ii,jj)=OutMAP.times(end);
        GNM_CPU(ii,jj)=OutGNM.times(end);

        RK_iter(ii,jj)=OutRK.iter;
        PIM_iter(ii,jj)=OutPIM.iter;
        MAP_iter(ii,jj)=OutMAP.iter;
        GNM_iter(ii,jj)=OutGNM.iter;
        %%%
        fprintf('Average time =%d\n',ii);
    end
    fprintf('m=%d \n',m);
    fprintf('CPU: PIM, GNM, MAP, RK \n');
    fprintf(' %4.2f&  %4.2f &  %4.2f &  %4.2f \n \n', mean(PIM_CPU(:,jj)),mean(GNM_CPU(:,jj)),...
        mean(MAP_CPU(:,jj)),mean(RK_CPU(:,jj)));
end


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
xlable=m1*[1:Numbn];

%%%
y1=RK_CPU';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%%
y2=PIM_CPU';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);


%%%
y3=MAP_CPU';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);



%%%
y4=GNM_CPU';
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
    'LineStyle', '-', 'Marker', 'o','DisplayName', 'RIMCS');
p3=semilogy( xlable, median(y3'), 'red', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker', 's','DisplayName', 'RIMGS');
p4=semilogy( xlable, median(y4'), 'green', 'LineWidth', 1,...
    'LineStyle', '-', 'Marker', '^','DisplayName', 'RIMSRHT');

p1=semilogy( xlable, median(y1'), 'magenta', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 'x', 'DisplayName', 'RABK');
set(gca, 'YScale', 'log')
ylabel('RSE','Interpreter', 'latex')
xlabel('Number of iterations','Interpreter', 'latex')
legend([ p2 p3 p4 p1],{'PIM','MAP','GNM','RABK'},'Interpreter', 'latex','location', 'best')
txt=title(['$\kappa_A=$ ',num2str(kappaA),', $\kappa_B=$ ',num2str(kappaB)]);
set(txt, 'Interpreter', 'latex');




