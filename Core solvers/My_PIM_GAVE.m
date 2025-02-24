function [x,Out]=My_PIM_GAVE(A,B,b,opts)

% Picard iteration method for solving GAVE
%              Ax-B|x|=b
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
% [1] Mangasarian, O.L.: A generalized Newton method for absolute value equation.
% Optim. Lett. 3, 101–108 (2008)
%
% Coded by Jiaxin Xie, Beihang University, xiejx@buaa.edu.cn

tic

%% setting some parameters
[m,n]=size(A);

%if m~=n
%    fprintf('the coefficient matrices must be square');
%    return;
%end

flag=exist('opts');

%%%% setting the max iteration
if (flag && isfield(opts,'Max_iter'))
    Max_iter=opts.Max_iter;
else
    Max_iter=1000;
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
    error1=norm(A*x-B*abs(x)-b)^2/normb;
end

RSE(1)=error1;
%[L,U,P] = lu(A);

%% executing the AmRABK method
stopc=0;
iter=0;
times(1)=toc;
while ~stopc
    tic
    iter=iter+1;

    %%
    %xold=x;
    %y = L\(P*(B*abs(x)+b));
    %x = U\y;
    x=A\(B*abs(x)+b);

    %% stopping rule
    if strategy
        error1=norm(x-xstar)^2/normxstar;
        RSE(iter+1)=error1;% RSE
        if error1<TOL  || iter>=Max_iter
            stopc=1;
        end
    else
        %%%% Note that we do not us this stopping rule during our test
        error1=norm(A*x-B*abs(x)-b)^2/normb;
        RSE(iter+1)=error1;
        if  error1<TOL || iter>=Max_iter
            stopc=1;
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

