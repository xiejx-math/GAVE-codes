function [x,Out]=My_RIMGS_GAVE(A,B,b,alpha,ell,opts)

% Randomized iterative method with Guassian Sketching  solving GAVE
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
% if (flag && isfield(opts,'permS'))
%     S=opts.permS;
%     A=A(S,:);
%     B=B(S,:);
%     b=b(S);
% else
%     S=randperm(m);
%     A=A(S,:);
%     B=B(S,:);
%     b=b(S);
% end

%% setting the probability
% if (flag && isfield(opts,'probset'))
%     probset=opts.probset;
% else
%     probset=0;
% end
% 
% if probset
%     Aarrs=opts.Aarrs;
%     Barrs=opts.Barrs;
%     barrs=opts.barrs;
%     cumsumpro=opts.cumsumpro;
% else
% 
%     normAfro=norm(A,'fro')^2;
%     tau=floor(m/ell);
%     blockAnormfro=zeros(tau,1);
%     %prob=zeros(tau,1);
%     for i=1:tau
%         if i==tau
%             ps=((i-1)*ell+1):1:m;
%         else
%             ps=((i-1)*ell+1):1:(i*ell);
%         end
%         Aps=A(ps,:);
%         blockAnormfro(i)=norm(A(ps,:),'fro')^2;
%         Aarrs{i}=Aps;
%         Barrs{i}=B(ps,:);
%         barrs{i}=b(ps);
%     end
%     prob=blockAnormfro/normAfro;
%     cumsumpro=cumsum(prob);
% end

%% executing the AmRABK method
stopc=0;
iter=0;
times(1)=toc;
while ~stopc
    tic
    iter=iter+1;

    %%

    %indexR=randperm(m,ell);
    %AindexR=A(indexR,:);
    %BindexR=B(indexR,:);
    %bindexR=b(indexR);

    [AindexR, BindexR, bindexR] = My_Gaussian_sketch(A, B, b, ell);
  
    Ax_Bx_b=AindexR*x-BindexR*abs(x)-bindexR;
    normAindexR=norm(AindexR);
    x=x-alpha*AindexR'*Ax_Bx_b/normAindexR^2;


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

