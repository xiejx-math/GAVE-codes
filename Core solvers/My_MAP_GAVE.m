function [x,Out]=My_MAP_GAVE(A,B,b,opts)

% methods of alternating projections for solving GAVE
%              Ax-B|x|=b
%
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
% Based on the paper:
% [1] J.H. Alcantara, J.-S. Chen, and M. K. Tam. Method of alternating projections for
% the general absolute value equation. Journal of Fixed Point Theory and Applications, 25(1):39, 2023
%
% Coded by Jiaxin Xie, Beihang University, xiejx@buaa.edu.cn

tic

%% setting some parameters
[m,n]=size(A);

flag=exist('opts');

%%%% setting the max iteration
if (flag && isfield(opts,'Max_iter'))
    Max_iter=opts.Max_iter;
else
    Max_iter=10000;
end

%%%% setting the tolerance
if (flag && isfield(opts,'TOL'))
    TOL=opts.TOL;
else
    TOL=10^-12;
end

%%%% dense or sparse data
if (flag && isfield(opts,'sparse'))
    sparsedata=opts.sparse;
else
    sparsedata=0;
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

w=[x;x];
%%
T=[A-B -A-B];
sqrt2c=sqrt(2)*b;

if sparsedata
    [U, S, V] = svd(T, 'econ');
else
    pinvT=pinv(T);
end
%S=2*(A*A'+B*B');
%dS=decomposition(S,'chol');

%norm(T*(pinvT*b)-b)

%% executing the AP method
stopc=0;
iter=0;
times(1)=toc;
while ~stopc
    tic
    iter=iter+1;

    %%
    w=projectionC2(w);
    if sparsedata
        %w=w-T\(T*w-sqrt2c);
        %norm(lsqminnorm(T,T*w-sqrt2c)-R(1:n,1:n)\(Q'*(T*w-sqrt2c)))
        %xxx(p,:) = R\(Q\(T*w-sqrt2c));
        
        w=w-lsqminnorm(T,T*w-sqrt2c);
       % w=w-R\(Q'*(T*w-sqrt2c));
       norm(www-w)
    else
        w=w-pinvT*(T*w-sqrt2c);
    end
    x=1/sqrt(2)*(w(1:n)-w(n+1:2*n));

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

