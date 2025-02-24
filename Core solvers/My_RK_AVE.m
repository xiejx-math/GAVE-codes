function [x,Out]=My_RK_AVE(A,b,alpha,opts)

% Randomized Kaczmarz for solving GAVE
%              Ax-|x|=b
% we use a simple partitioning strategy for chosing the sampling matrix
%
%Input: the coefficent matrices A and B, the vector b and opts
%opts.initial: the initial vector x^0
%opts.TOL: the stopping rule
%.....
%
%Output: the approximate solution x and Out
% Out.error: the relative iterative error \|x^k-x^*\|^2/\|x^k\|^2
% Out.iter: the total number of iteration
% ....
%
% Based on the manuscript:
% [1] Jiaxin Xie, Houduo Qi, and Deren Han. Randomized iterative methods for generalized absolute value equations: Solvability and error bounds,
%  arXiv:\
%
% Coded by Jiaxin Xie, Beihang University, xiejx@buaa.edu.cn

tic

%% setting some parameters
[m,n]=size(A);
if m~=n
    fprintf('the coefficient matrices must be square');
    return;
end

flag=exist('opts');

%%%% setting the max iteration
if (flag && isfield(opts,'Max_iter'))
    Max_iter=opts.Max_iter;
else
    Max_iter=2000000;
end

%%%% setting the tolerance
if (flag && isfield(opts,'TOL'))
    TOL=opts.TOL;
else
    TOL=10^-12;
end

%%%% setting the initial point
if (flag && isfield(opts,'initial'))
    initialx=opts.initial;
else
    initialx=zeros(n,1);
end
x=initialx;


%%%% determining what to use as the stopping rule
if (flag && isfield(opts,'xstar'))
    xstar=opts.xstar;
    if m>=n
        normxstar=norm(xstar)^2;
        error1=norm(xstar-x)^2/normxstar;
        strategy=1;
    else
        strategy=0;
    end
else
    strategy=0;
end

if (flag && isfield(opts,'strategy'))
    strategy=opts.strategy;
    normxstar=norm(xstar)^2;
end

if ~strategy
    normb=norm(b)^2+1;
    error1=norm(A*x-abs(x)-b)^2/normb;
end

RSE(1)=error1;


%% setting the probability
Anrm2=sum(A.^2,2);
pro1=Anrm2/norm(A,'fro')^2;
cumsumpro=cumsum(pro1);
mm=min(100,m);

l1=sum(cumsumpro<rand(1,mm),1)+1;
iter1=0;
iter2=0;

%% executing the AmRABK method
stopc=0;
iter=0;
times(1)=toc;
while ~stopc
    tic
    iter=iter+1;

    %%
    iter1=iter1+1;
    if iter1>mm
        iter1=1;
        l1=sum(cumsumpro<rand(1,mm),1)+1;
    end
    %%

    l=l1(iter1);

    %l=sum(cumsumpro<rand)+1;
    AindexR=A(l,:);
    AindexRt=AindexR';
    %BindexR=B(l,:);
    bindexR=b(l);
    Ax_Bx_b=AindexR*x-abs(x(l))-bindexR;
    %x=x-alpha*AindexR'*Ax_Bx_b/Anrm2(l);
    alphaAnormAb=alpha*Ax_Bx_b/Anrm2(l);
    %AtindexRAxBxb=AindexR';

    
    x=x-alphaAnormAb*AindexRt;

    %% stopping rule
    if strategy
        error1=norm(x-xstar)^2/normxstar;
        RSE(iter+1)=error1;% RSE
        if error1<TOL  || iter>=Max_iter
            stopc=1;
        end
    else
        iter2=iter2+1;
        if iter2>mm
            %%%% Note that we do not us this stopping rule during our test
            iter2=1;
            error1=norm(A*x-abs(x)-b)^2/normb;
            RSE(iter+1)=error1;
            if  error1<TOL || iter>=Max_iter
                stopc=1;
            end
        end
    end
    %% store the CPU time
    times(iter+1)=times(iter)+toc;
end
%% setting Output
Out.error=RSE;
Out.iter=iter;
Out.times=times;
end

