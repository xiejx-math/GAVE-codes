% This Matlab file is used to compare PIM, GNM, MAP,
% and RABK with square matrices.
% This file can reproduce Figure 2 in the manuscript

clear
close all

m1=1000;
dense=1;
run_time=20;
Numbn=3;

CPUdata=zeros(4,Numbn);


for jj=1:Numbn
    m=m1*jj;
    n=m;
    r=min(m,n);
    RK_CPU=zeros(run_time,1);
    PIM_CPU=zeros(run_time,1);
    MAP_CPU=zeros(run_time,1);
    GNM_CPU=zeros(run_time,1);

    RK_iter=zeros(run_time,1);
    PIM_iter=zeros(run_time,1);
    MAP_iter=zeros(run_time,1);
    GNM_iter=zeros(run_time,1);

    for ii=1:run_time
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
        RK_CPU(ii)=OutRK.times(end);
        PIM_CPU(ii)=OutPIM.times(end);
        MAP_CPU(ii)=OutMAP.times(end);
        GNM_CPU(ii)=OutGNM.times(end);

        RK_iter(ii)=OutRK.iter;
        PIM_iter(ii)=OutPIM.iter;
        MAP_iter(ii)=OutMAP.iter;
        GNM_iter(ii)=OutGNM.iter;
        %%
        fprintf('Average time =%d\n',ii)
    end

    fprintf('m=%d \n',m)
    fprintf('Iter and CPU: PIM, GNM, MAP, RK \n')
    fprintf(' %4.2f&  %4.2f &  %4.2f &  %4.2f & %4.2f&  %4.2f &  %4.2f &  %4.2f\n',mean(PIM_iter),...
        mean(PIM_CPU),mean(GNM_iter),mean(GNM_CPU),mean(MAP_iter),mean(MAP_CPU),...
        mean(RK_iter),mean(RK_CPU))
    CPUdata(1,jj)=mean(PIM_CPU);
    CPUdata(2,jj)=mean(GNM_CPU);
    CPUdata(3,jj)=mean(MAP_CPU);
    CPUdata(4,jj)=mean(RK_CPU);

end



x_label=m1*[1:Numbn];

figure
plot(x_label,CPUdata(1,:),'red', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 'o')

hold on
plot(x_label,CPUdata(2,:),'blue', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 's')
plot(x_label,CPUdata(3,:),'green', 'LineWidth', 1,...
    'LineStyle', '-','Marker', '^')
plot(x_label,CPUdata(4,:),'magenta', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 'x')
ylabel('CPU')
xlabel('The values of $n$','Interpreter', 'latex')
legend('PIM','GNM','MAP','RABK','location', 'best')

figure
semilogy(x_label,CPUdata(1,:),'red', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 'o')

hold on
semilogy(x_label,CPUdata(2,:),'blue', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 's')
semilogy(x_label,CPUdata(3,:),'green', 'LineWidth', 1,...
    'LineStyle', '-','Marker', '^')
semilogy(x_label,CPUdata(4,:),'magenta', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 'x')
ylabel('CPU')
xlabel('The values of $n$','Interpreter', 'latex')
legend('PIM','GNM','MAP','RABK','location', 'best')



