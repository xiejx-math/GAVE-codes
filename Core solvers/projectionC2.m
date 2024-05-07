function w=projectionC2(w)

n2=length(w);

u=w(1:n2/2);
v=w((n2/2+1):n2);
%%
I=find(u<=v);
u(I)=0;
v(I)=max(v(I),0);
%%
J=find(u>v);%setdiff([1:n2/2],I);
v(J)=0;
u(J)=max(u(J),0);

w=[u;v];

end