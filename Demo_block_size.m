% This Matlab file is used to explore the influence of the block size
% This file can reproduce Figure 2 in the manuscript

clear
close all

ss=9;
m=2^ss;
kappaA=2;
kappaB=10;



%%
setvaluep=[0:1:log2(m)];
sizeP=length(setvaluep);
k_time=sizeP;
run_time=20;


%% n=256
n=2^(ss-2);
r=min(m,n);

RABKtimes256=zeros(run_time,k_time);
RABKiter256=zeros(run_time,k_time);
for ii=1:k_time
    ell=2^(setvaluep(ii));
    for jj=1:run_time
        %%
        %B=randn(m,n);
        %B=eye(n);
        %[U,~]=qr(randn(m,r),0);
        %[V,~]=qr(randn(n,r),0);
        %sigmaA=norm(B)*eye(r)+diag(rand(r,1));
        %A=U*sigmaA*V';
        %%%
        [A,B]=getAB(m,n,2,kappaA,1,kappaB);

        %%%
        x=randn(n,1);
        b=A*x-B*abs(x);

        %%
        opts.xstar=x;

        alpha=1;
        [xRABK,OutRABK]=My_RABK_GAVE(A,B,b,alpha,ell,opts);

        %%
        RABKtimes256(jj,ii)=OutRABK.times(end);
        RABKiter256(jj,ii)=OutRABK.iter;

    end
    fprintf('ell=%d\n',ell)
end



%% n=512
n=2^(ss-1);
r=min(m,n);

RABKtimes512=zeros(run_time,k_time);
RABKiter512=zeros(run_time,k_time);
for ii=1:k_time
    ell=2^(setvaluep(ii));
    for jj=1:run_time
        %B=randn(m,n);
        %B=eye(n);
        %[U,~]=qr(randn(m,r),0);
        %[V,~]=qr(randn(n,r),0);
        %sigmaA=norm(B)*eye(r)+diag(rand(r,1));
        %A=U*sigmaA*V';

        %%%
        [A,B]=getAB(m,n,2,kappaA,1,kappaB);
        x=randn(n,1);
        b=A*x-B*abs(x);

        %%
        opts.xstar=x;

        alpha=1;
        [xRABK,OutRABK]=My_RABK_GAVE(A,B,b,alpha,ell,opts);

        %%
        RABKtimes512(jj,ii)=OutRABK.times(end);
        RABKiter512(jj,ii)=OutRABK.iter;

    end
    fprintf('ell=%d\n',ell)
end

%% n=1024
n=2^ss;
r=min(m,n);

RABKtimes1024=zeros(run_time,k_time);
RABKiter1024=zeros(run_time,k_time);
for ii=1:k_time
    ell=2^(setvaluep(ii));
    for jj=1:run_time
        %B=randn(m,n);
        %B=eye(n);
        %[U,~]=qr(randn(m,r),0);
        %[V,~]=qr(randn(n,r),0);
        %sigmaA=norm(B)*eye(r)+diag(rand(r,1));
        %A=U*sigmaA*V';

        %%%
        [A,B]=getAB(m,n,2,kappaA,1,kappaB);

        %%%
        x=randn(n,1);
        b=A*x-B*abs(x);

        %%
        opts.xstar=x;

        alpha=1;
        [xRABK,OutRABK]=My_RABK_GAVE(A,B,b,alpha,ell,opts);

        %%
        RABKtimes1024(jj,ii)=OutRABK.times(end);
        RABKiter1024(jj,ii)=OutRABK.iter;

    end
    fprintf('ell=%d\n',ell)
end

%%

num_iter_array=setvaluep';
P_matrix=zeros(run_time,k_time);
for kk=1:run_time
    P_matrix(kk,:)=2.^(num_iter_array');
end

RABK_value256=RABKiter256.*P_matrix/m;
RABK_value512=RABKiter512.*P_matrix/m;
RABK_value1024=RABKiter1024.*P_matrix/m;

%%
xlable=setvaluep;

%%
y1=RABKtimes256';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%
y2=RABKtimes512';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%
y3=RABKtimes1024';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);

%%

figure
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
p1=semilogy( xlable, median(y1'), 'red', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)

p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
            'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
ylabel('CPU')
xlabel('$\log_2(p)$','Interpreter', 'latex')
%legend('RABK','location', 'best')
legend([p1 p2 p3],{'RABK $n=128$','RABK $n=256$','RABK $n=512$'},'Interpreter', 'latex','location', 'best')
%txt=title(['{\tt randn}',',$m=$ ',num2str(m)]);
txt=title(['$m=$ ',num2str(m),', $\kappa_A=$ ',num2str(kappaA),', $\kappa_B=$ ',num2str(kappaB)]);
set(txt, 'Interpreter', 'latex');

%%
xlable=setvaluep;

%%
y1=RABK_value256';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%
y2=RABK_value512';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%
y3=RABK_value1024';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);

figure
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .2)
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .2)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .2)
p1=semilogy( xlable, median(y1'), 'red', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
            'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
ylabel('$k\cdot \frac{p}{m}$','Interpreter', 'latex')
xlabel('$\log_2(p)$','Interpreter', 'latex')
legend([p1 p2 p3],{'RABK $n=128$','RABK $n=256$','RABK $n=512$'},'Interpreter', 'latex','location', 'best')
%txt=title(['{\tt randn}',',$m=$ ',num2str(m)]);
txt=title(['$m=$ ',num2str(m),', $\kappa_A=$ ',num2str(kappaA),', $\kappa_B=$ ',num2str(kappaB)]);
set(txt, 'Interpreter', 'latex');
