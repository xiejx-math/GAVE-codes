function [x,Out]=My_RABK_GAVE(A,B,b,alpha,ell,opts)

% Randomized average block Kaczmarz for solving GAVE
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
% [1] Jiaxin Xie, Houduo Qi, and Deren Han. Randomized iterative methods for generalized absolute value equations: Solvability and error bounds,
%  arXiv:\
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
    Max_iter=2000000;
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

%% a uniform random permutation for both A and b
if (flag && isfield(opts,'permS'))
    S=opts.permS;
    A=A(S,:);
    B=B(S,:);
    b=b(S);
else
    S=randperm(m);
    A=A(S,:);
    B=B(S,:);
    b=b(S);
end

%% setting the probability
if (flag && isfield(opts,'probset'))
    probset=opts.probset;
else
    probset=0;
end

if probset
    Aarrs=opts.Aarrs;
    Barrs=opts.Barrs;
    barrs=opts.barrs;
    cumsumpro=opts.cumsumpro;
else

    %normAfro=norm(A,'fro')^2;
    tau=floor(m/ell);
    blockAnormfro=zeros(tau,1);
    %prob=zeros(tau,1);
    for i=1:tau
        if i==tau
            ps=((i-1)*ell+1):1:m;
        else
            ps=((i-1)*ell+1):1:(i*ell);
        end
        Aps=A(ps,:);
        if sparsedata
            blockAnormfro(i)=normest(A(ps,:),1.0e-4)^2;
        else
            blockAnormfro(i)=norm(A(ps,:))^2;
        end
        Aarrs{i}=Aps;
        Barrs{i}=B(ps,:);
        barrs{i}=b(ps);
    end
    prob=blockAnormfro/sum(blockAnormfro);
    cumsumpro=cumsum(prob);
end



%% executing the AmRABK method
stopc=0;
iter=0;
times(1)=toc;
while ~stopc
    tic
    iter=iter+1;

    %%
    l=sum(cumsumpro<rand)+1;
    AindexR=Aarrs{l};
    BindexR=Barrs{l};
    bindexR=barrs{l};
    Ax_Bx_b=AindexR*x-BindexR*abs(x)-bindexR;
    x=x-alpha*AindexR'*Ax_Bx_b/blockAnormfro(l);

%     %% Randomly select a sampling matrix S_k until $S_k(A*x_k − b) \neq 0$
%     stopsampling=0;
%     while ~stopsampling
%         l=sum(cumsumpro<rand)+1;
%         AindexR=Aarrs{l};
%         bindexR=barrs{l};
%         Axb=AindexR*x-bindexR;
%         normAxb=norm(Axb)^2;
% 
%         % If normAxb is greater than 10^(-16), we consider it to be equal to 0
%         if normAxb>TOL1
%             stopsampling=1;
%             dk=AindexR'*Axb;
%             norm_dk=norm(dk)^2;
%         end
%     end
% 
%     %% compute the x
%     if iter==1
%         alpha=normAxb/norm_dk;
%         xold=x;
%         %%%% update x for k=1
%         x=x-alpha*dk;
%     else
%         x_xoold=x-xold;
%         %normAxb=norm(Axb)^2;
%         norm_x_xoold=norm(x_xoold)^2;
%         dk_x_xoold=dk'*x_xoold;
%         denomfm=norm_dk*norm_x_xoold-dk_x_xoold^2;
%         alpha=(normAxb*norm_x_xoold)/denomfm;
%         beta=(dk_x_xoold*normAxb)/denomfm;
%         %%%% update x for k \geq 2
%         xoold=xold;
%         xold=x;
%         x=x-alpha*dk+beta*(x-xoold);
%     end

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

