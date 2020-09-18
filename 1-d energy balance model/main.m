clear;
clc;

global c dab1 Ch Gamma theta S0
Ch=46.8;
Gamma=0.61;
theta=5.67e-8;
S0=1368;
c=1;

fxy=@(x)(S0*(0.5+0.2*tanh((x-265)/10))/4-Gamma*theta*x^4)/Ch;
xnode1=fsolve(fxy,225);
xnode2=fsolve(fxy,300);
xsad=fsolve(fxy,260);
xnode=xnode1;

% syms x1f p1f c dab1 Ch Gamma theta S0
% f1=(S0*(0.5+0.2*tanh((x1f-265)/10))/4-Gamma*theta*x1f^4)/Ch;
% f1x1=gradient(f1,x1f);
% H=c^2*p1f^2+2*(f1-dab1)*p1f-f1x1+dab1^2/c^2;
% Hx1=gradient(H,x1f);
% Hp1=gradient(H,p1f);

alpha1=0.5;
beta1=0;
dab1=alpha1*beta1/gamma(2-alpha1)/cos(pi*alpha1/2);

%%% ³õÖµ
m=1e3;
Npath=1e3;
T=linspace(0,1,Npath+1);
h=T(2)-T(1);
I=ones(m,1);
x1=xnode(1)+0.1*I*T;
p1=zeros(m,Npath+1);
lambdamax=100;
lambda1=2*lambdamax*rand(m,1)-lambdamax;
% lam=5;
% lambda1=lam*randn(m,1);
% lambda2=lam*randn(m,1);
p1(:,end)=-lambda1;
x01=x1;
p01=p1;
delta=0.0001;
xT1=zeros(m,1);
pos=1:m;
Index=0;

%%% µü´ú
for i=1:10000000
    for k=1:Npath
        x1f=x01(:,Npath+2-k);
        p1f=p01(:,Npath+2-k);
        K1p1=h*(- ((S0*tanh(x1f/10 - 53/2).*(tanh(x1f/10 - 53/2).^2/10 - 1/10))/100 - 12*Gamma*theta*x1f.^2)/Ch -...
            (2*p1f.*((S0*(tanh(x1f/10 - 53/2).^2/50 - 1/50))/4 + 4*Gamma*theta*x1f.^3))/Ch);
        
        x1f=x01(:,Npath+1-k);
        p1f=p01(:,Npath+2-k)+K1p1;
        K2p1=h*(- ((S0*tanh(x1f/10 - 53/2).*(tanh(x1f/10 - 53/2).^2/10 - 1/10))/100 - 12*Gamma*theta*x1f.^2)/Ch -...
            (2*p1f.*((S0*(tanh(x1f/10 - 53/2).^2/50 - 1/50))/4 + 4*Gamma*theta*x1f.^3))/Ch);
        
        p1(:,Npath+1-k)=p01(:,Npath+2-k)+1/2*(K1p1+K2p1);
    end
    
    for k=1:Npath
        x1f=x01(:,k);
        p1f=p01(:,k);
        K1x1=h*(2*p1f*c^2 - 2*dab1 + ((S0*(tanh(x1f/10 - 53/2)/5 + 1/2))/2 - 2*Gamma*theta*x1f.^4)/Ch);
        
        x1f=x01(:,k)+K1x1;
        p1f=p01(:,k+1);
        K2x1=h*(2*p1f*c^2 - 2*dab1 + ((S0*(tanh(x1f/10 - 53/2)/5 + 1/2))/2 - 2*Gamma*theta*x1f.^4)/Ch);
        
        x1(:,k+1)=x01(:,k)+1/2*(K1x1+K2x1);
    end
    
    K=length(x1(:,1));
    I=[];
    for k=1:K
        if (max(abs(x1(k,:)-x01(k,:)))<delta)&&(i>1000)
            xT1(pos(k))=x1(k,end);
            I=[I k];
        end
    end
    if ~isempty(I)
        pos(I)=[];
        x1(I,:)=[];
        p1(I,:)=[];
    end
    
    if isempty(x1)
        break;
    end
    x01=x1;
    p01=p1;
end

xT0=xT1;
xT1=2/(xnode2-xnode1)*(xT1-(xnode1+xnode2)/2);
lambdascale=100;
lambda0=lambda1;
lambda1=lambda1/lambdascale;

figure;
plot(xT0,'*');

figure;
plot(lambda0,'.');

figure;
plot(xT0,lambda0,'*');

% figure;
% plot([0 1000],[xnode1 xnode1]);
% hold on
% plot([0 1000],[xnode2 xnode2]);
% hold off

%%% Neural Network
eta=0.01;
A0=xT1';
Lambda=lambda1';
n1=20;
n2=20;
n3=20;
nf=1;
DELTAw=0.1;
DELTAb=0.01;
W1=DELTAw*randn(n1,1);
b1=DELTAb*randn(n1,1);
W2=DELTAw*randn(n2,n1);
b2=DELTAb*randn(n2,1);
W3=DELTAw*randn(n3,n2);
b3=DELTAb*randn(n3,1);
Wf=DELTAw*randn(nf,n3);
bf=DELTAb*randn(nf,1);
% Wf=DELTAw*randn(nf,n1);
% bf=DELTAb*randn(nf,1);
M=1/m*ones(m,1);
I=ones(1,m);
Jcost=[];

for i=1:200000
    Z1=W1*A0+b1*I;
    A1=Z1;
    I0=A1<0;
    A1(I0)=0;
    
    Z2=W2*A1+b2*I;
    A2=Z2;
    I0=A2<0;
    A2(I0)=0;
    
    Z3=W3*A2+b3*I;
    A3=Z3;
    I0=A3<0;
    A3(I0)=0;
    
    Zf=Wf*A3+bf*I;
    Af=Zf;
    
%     Zf=Wf*A1+bf*I;
%     Af=Zf;
    
    Jcost=[Jcost sum((Af-Lambda).^2)/2/m];
    
    Ldaf=Af-Lambda;
    Ldzf=Ldaf;
%     Jdwf=Ldzf*A1'/m;
%     Jdbf=Ldzf*M;
    Jdwf=Ldzf*A3'/m;
    Jdbf=Ldzf*M;
    
    Lda3=Wf'*Ldzf;
    Ldz3=Lda3.*(Z3>=0);
    Jdw3=Ldz3*A2'/m;
    Jdb3=Ldz3*M;
    
    Lda2=W3'*Ldz3;
    Ldz2=Lda2.*(Z2>=0);
    Jdw2=Ldz2*A1'/m;
    Jdb2=Ldz2*M;
    
    Lda1=W2'*Ldz2;
    Ldz1=Lda1.*(Z1>=0);
    Jdw1=Ldz1*A0'/m;
    Jdb1=Ldz1*M;
    
%     Lda1=Wf'*Ldzf;
%     Ldz1=Lda1.*(Z1>=0);
%     Jdw1=Ldz1*A0'/m;
%     Jdb1=Ldz1*M;
    
    W1=W1-eta*Jdw1;
    b1=b1-eta*Jdb1;
    W2=W2-eta*Jdw2;
    b2=b2-eta*Jdb2;
    W3=W3-eta*Jdw3;
    b3=b3-eta*Jdb3;
    Wf=Wf-eta*Jdwf;
    bf=bf-eta*Jdbf;
    
    if max(max(abs(Jdw1)))+max(abs(Jdb1))+max(max(abs(Jdw2)))+max(abs(Jdb2))+max(max(abs(Jdw3)))+max(abs(Jdb3))+max(max(abs(Jdwf)))+max(abs(Jdbf))<0.0001
%     if max(max(abs(Jdw1)))+max(abs(Jdb1))+max(max(abs(Jdwf)))+max(abs(Jdbf))<0.01
        break;
    end
end

figure;
plot(Jcost);

%%% Test
xend=1;
xf=xend;
z1=W1*xf+b1;
a1=z1;
I0=a1<0;
a1(I0)=0;

z2=W2*a1+b2;
a2=z2;
I0=a2<0;
a2(I0)=0;

z3=W3*a2+b3;
a3=z3;
I0=a3<0;
a3(I0)=0;

zf=Wf*a3+bf;
af=zf;

af=af*lambdascale;

% zf=Wf*a1+bf;
% af=zf;

Npath=1e3;
x1=xnode(1)+0.1*T;
p1=zeros(1,Npath+1);
p1(end)=-af;
x01=x1;
p01=p1;
delta=0.0001;

for i=1:100000
    for k=1:Npath
        x1f=x01(:,Npath+2-k);
        p1f=p01(:,Npath+2-k);
        K1p1=h*(- ((S0*tanh(x1f/10 - 53/2)*(tanh(x1f/10 - 53/2)^2/10 - 1/10))/100 - 12*Gamma*theta*x1f^2)/Ch -...
            (2*p1f*((S0*(tanh(x1f/10 - 53/2)^2/50 - 1/50))/4 + 4*Gamma*theta*x1f^3))/Ch);
        
        x1f=x01(:,Npath+1-k);
        p1f=p01(:,Npath+2-k)+K1p1;
        K2p1=h*(- ((S0*tanh(x1f/10 - 53/2)*(tanh(x1f/10 - 53/2)^2/10 - 1/10))/100 - 12*Gamma*theta*x1f^2)/Ch -...
            (2*p1f*((S0*(tanh(x1f/10 - 53/2)^2/50 - 1/50))/4 + 4*Gamma*theta*x1f^3))/Ch);
        
        p1(:,Npath+1-k)=p01(:,Npath+2-k)+1/2*(K1p1+K2p1);
    end
    
    for k=1:Npath
        x1f=x01(:,k);
        p1f=p01(:,k);
        K1x1=h*(2*p1f*c^2 - 2*dab1 + ((S0*(tanh(x1f/10 - 53/2)/5 + 1/2))/2 - 2*Gamma*theta*x1f^4)/Ch);
        
        x1f=x01(:,k)+K1x1;
        p1f=p01(:,k+1);
        K2x1=h*(2*p1f*c^2 - 2*dab1 + ((S0*(tanh(x1f/10 - 53/2)/5 + 1/2))/2 - 2*Gamma*theta*x1f^4)/Ch);
        
        x1(:,k+1)=x01(:,k)+1/2*(K1x1+K2x1);
    end
    
    if (max(abs(x1-x01))<delta)&&(i>1000)
        break;
    end
    x01=x1;
    p01=p1;
end

figure;
plot(T,x1);

figure;
plot(T,p1);

%%% Validation
X=zeros(2,Npath+1);
X(:,1)=[x1(1);p1(1)];
t0=0;
for i=1:Npath
    X(:,i+1)=rk4(t0,h,X(:,i));
end

figure;
plot(T,X(1,:));

figure;
plot(T,X(2,:));


