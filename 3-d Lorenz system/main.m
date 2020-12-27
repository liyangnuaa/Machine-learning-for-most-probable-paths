clear;
clc;

global dab1 dab2 dab3 sigma beta rou kapa miu
sigma=1;
beta=8/3;
rou=0.5;
kapa=0.5;
miu=0;
xnode=[0;0;0];

% syms x1f x2f x3f p1f p2f p3f dab1 dab2 dab3 sigma beta rou kapa miu
% f1=sigma*(x2f-x1f);
% f2=rou*x1f-x2f-x1f*x3f;
% f3=-beta*x3f+x1f*x2f+(kapa-0.5)*miu*x3f;
% f1x1=gradient(f1,x1f);
% f2x2=gradient(f2,x2f);
% f3x3=gradient(f3,x3f);
% H=1/4*(p1f^2+p2f^2+(1+miu*x3f^2)*p3f^2)+(f1-dab1)*p1f+(f2-dab2)*p2f+(f3-dab3)*p3f-f1x1-f2x2-f3x3+miu*x3f*f3/(1+miu*x3f^2)+dab1^2+dab2^2+dab3^2/(1+miu*x3f^2);
% Hx1=gradient(H,x1f);
% Hx2=gradient(H,x2f);
% Hx3=gradient(H,x3f);
% Hp1=gradient(H,p1f);
% Hp2=gradient(H,p2f);
% Hp3=gradient(H,p3f);
%%%% vectorize(Hx1)

alpha1=0.5;
beta1=0;
dab1=alpha1*beta1/gamma(2-alpha1)/cos(pi*alpha1/2);
alpha2=0.5;
beta2=0;
dab2=alpha2*beta2/gamma(2-alpha2)/cos(pi*alpha2/2);
alpha3=0.5;
beta3=0;
dab3=alpha3*beta3/gamma(2-alpha3)/cos(pi*alpha3/2);

%%% 初值
m=1e3;
Npath=1e2;
T=linspace(0,1,Npath+1);
h=T(2)-T(1);
I=ones(m,1);
x1=xnode(1)+0.1*I*T;
x2=xnode(2)+0.1*I*T;
x3=xnode(2)+0.1*I*T;
p1=zeros(m,Npath+1);
p2=zeros(m,Npath+1);
p3=zeros(m,Npath+1);
lambdamax=10;
lambda1=2*lambdamax*rand(m,1)-lambdamax;
lambda2=2*lambdamax*rand(m,1)-lambdamax;
lambda3=2*lambdamax*rand(m,1)-lambdamax;
% lam=5;
% lambda1=lam*randn(m,1);
% lambda2=lam*randn(m,1);
p1(:,end)=lambda1;
p2(:,end)=lambda2;
p3(:,end)=lambda3;
x01=x1;
x02=x2;
x03=x3;
p01=p1;
p02=p2;
p03=p3;
delta=0.0001;
xT1=zeros(m,1);
xT2=zeros(m,1);
xT3=zeros(m,1);
pos=1:m;
I0=[];
Index=0;


% m00=5;
% x1=x1(1:m00,:);
% x2=x2(1:m00,:);
% x3=x3(1:m00,:);
% p1=p1(1:m00,:);
% p2=p2(1:m00,:);
% p3=p3(1:m00,:);
% x01=x1;
% x02=x2;
% x03=x3;
% p01=p1;
% p02=p2;
% p03=p3;
% for i=1:10000000
%     for k=1:Npath
%         x1f=x01(:,Npath+2-k);
%         x2f=x02(:,Npath+2-k);
%         x3f=x03(:,Npath+2-k);
%         p1f=p01(:,Npath+2-k);
%         p2f=p02(:,Npath+2-k);
%         p3f=p03(:,Npath+2-k);
%         K1p1=h*(2*p3f.*x2f - 2*p1f*sigma + p2f.*(2*rou - 2*x3f));
%         K1p2=h*(2*p1f*sigma - 2*p2f + 2*p3f.*x1f);
%         K1p3=h*(- 2*beta*p3f - 2*p2f.*x1f);
%         
%         x1f=x01(:,Npath+1-k);
%         x2f=x02(:,Npath+1-k);
%         x3f=x03(:,Npath+1-k);
%         p1f=p01(:,Npath+2-k)+K1p1;
%         p2f=p02(:,Npath+2-k)+K1p2;
%         p3f=p03(:,Npath+2-k)+K1p3;
%         K2p1=h*(2*p3f.*x2f - 2*p1f*sigma + p2f.*(2*rou - 2*x3f));
%         K2p2=h*(2*p1f*sigma - 2*p2f + 2*p3f.*x1f);
%         K2p3=h*(- 2*beta*p3f - 2*p2f.*x1f);
%         
%         p1(:,Npath+1-k)=p01(:,Npath+2-k)+1/2*(K1p1+K2p1);
%         p2(:,Npath+1-k)=p02(:,Npath+2-k)+1/2*(K1p2+K2p2);
%         p3(:,Npath+1-k)=p03(:,Npath+2-k)+1/2*(K1p3+K2p3);
%     end
%     
%     for k=1:Npath
%         x1f=x01(:,k);
%         x2f=x02(:,k);
%         x3f=x03(:,k);
%         p1f=p01(:,k);
%         p2f=p02(:,k);
%         p3f=p03(:,k);
%         K1x1=h*(2*p1f*c^2 - 2*dab1 - 2*sigma*(x1f - x2f));
%         K1x2=h*(2*p2f*c^2 - 2*dab2 - 2*x2f + 2*rou*x1f - 2*x1f.*x3f);
%         K1x3=h*(2*p3f*c^2 - 2*dab3 - 2*beta*x3f + 2*x1f.*x2f);
%         
%         x1f=x01(:,k)+K1x1;
%         x2f=x02(:,k)+K1x2;
%         x3f=x03(:,k)+K1x3;
%         p1f=p01(:,k+1);
%         p2f=p02(:,k+1);
%         p3f=p03(:,k+1);
%         K2x1=h*(2*p1f*c^2 - 2*dab1 - 2*sigma*(x1f - x2f));
%         K2x2=h*(2*p2f*c^2 - 2*dab2 - 2*x2f + 2*rou*x1f - 2*x1f.*x3f);
%         K2x3=h*(2*p3f*c^2 - 2*dab3 - 2*beta*x3f + 2*x1f.*x2f);
%         
%         x1(:,k+1)=x01(:,k)+1/2*(K1x1+K2x1);
%         x2(:,k+1)=x02(:,k)+1/2*(K1x2+K2x2);
%         x3(:,k+1)=x03(:,k)+1/2*(K1x3+K2x3);
%     end
%     
%     if (max(max(abs(x1-x01)+abs(x2-x02)+abs(x3-x03)))<delta)&&(i>1000)
%         break;
%     end
% 
%     x01=x1;
%     x02=x2;
%     x03=x3;
%     p01=p1;
%     p02=p2;
%     p03=p3;
% end
% figure;
% for i=1:m00
%     plot(T,x1(i,:));
%     hold on
% end
% hold off



%%% 迭代
for i=1:10000000
    for k=1:Npath
        x1f=x01(:,Npath+2-k);
        x2f=x02(:,Npath+2-k);
        x3f=x03(:,Npath+2-k);
        p1f=p01(:,Npath+2-k);
        p2f=p02(:,Npath+2-k);
        p3f=p03(:,Npath+2-k);
        K1p1=h*(p3f.*x2f - p1f.*sigma + p2f.*(rou - x3f) + (miu.*x2f.*x3f)./(miu.*x3f.^2 + 1));
        K1p2=h*(p1f.*sigma - p2f + p3f.*x1f + (miu.*x1f.*x3f)./(miu.*x3f.^2 + 1));
        K1p3=h*((miu.*p3f.^2.*x3f)./2 - p2f.*x1f - p3f.*(beta - miu.*(kapa - 1./2)) +...
            (miu.*(x1f.*x2f - beta.*x3f + miu.*x3f.*(kapa - 1./2)))./(miu.*x3f.^2 + 1) -...
            (miu.*x3f.*(beta - miu.*(kapa - 1./2)))./(miu.*x3f.^2 + 1) -...
            (2.*miu.^2.*x3f.^2.*(x1f.*x2f - beta.*x3f + miu.*x3f.*(kapa - 1./2)))./(miu.*x3f.^2 + 1).^2 -...
            (2.*dab3.^2.*miu.*x3f)./(miu.*x3f.^2 + 1).^2);
        
        x1f=x01(:,Npath+1-k);
        x2f=x02(:,Npath+1-k);
        x3f=x03(:,Npath+1-k);
        p1f=p01(:,Npath+2-k)+K1p1;
        p2f=p02(:,Npath+2-k)+K1p2;
        p3f=p03(:,Npath+2-k)+K1p3;
        K2p1=h*(p3f.*x2f - p1f.*sigma + p2f.*(rou - x3f) + (miu.*x2f.*x3f)./(miu.*x3f.^2 + 1));
        K2p2=h*(p1f.*sigma - p2f + p3f.*x1f + (miu.*x1f.*x3f)./(miu.*x3f.^2 + 1));
        K2p3=h*((miu.*p3f.^2.*x3f)./2 - p2f.*x1f - p3f.*(beta - miu.*(kapa - 1./2)) +...
            (miu.*(x1f.*x2f - beta.*x3f + miu.*x3f.*(kapa - 1./2)))./(miu.*x3f.^2 + 1) -...
            (miu.*x3f.*(beta - miu.*(kapa - 1./2)))./(miu.*x3f.^2 + 1) -...
            (2.*miu.^2.*x3f.^2.*(x1f.*x2f - beta.*x3f + miu.*x3f.*(kapa - 1./2)))./(miu.*x3f.^2 + 1).^2 -...
            (2.*dab3.^2.*miu.*x3f)./(miu.*x3f.^2 + 1).^2);
        
        p1(:,Npath+1-k)=p01(:,Npath+2-k)+1/2*(K1p1+K2p1);
        p2(:,Npath+1-k)=p02(:,Npath+2-k)+1/2*(K1p2+K2p2);
        p3(:,Npath+1-k)=p03(:,Npath+2-k)+1/2*(K1p3+K2p3);
    end
    
    for k=1:Npath
        x1f=x01(:,k);
        x2f=x02(:,k);
        x3f=x03(:,k);
        p1f=p01(:,k);
        p2f=p02(:,k);
        p3f=p03(:,k);
        K1x1=h*(p1f./2 - dab1 - sigma.*(x1f - x2f));
        K1x2=h*(p2f./2 - dab2 - x2f + rou.*x1f - x1f.*x3f);
        K1x3=h*(x1f.*x2f - beta.*x3f - dab3 + (p3f.*(miu.*x3f.^2 + 1))./2 + miu.*x3f.*(kapa - 1./2));
        
        x1f=x01(:,k)+K1x1;
        x2f=x02(:,k)+K1x2;
        x3f=x03(:,k)+K1x3;
        p1f=p01(:,k+1);
        p2f=p02(:,k+1);
        p3f=p03(:,k+1);
        K2x1=h*(p1f./2 - dab1 - sigma.*(x1f - x2f));
        K2x2=h*(p2f./2 - dab2 - x2f + rou.*x1f - x1f.*x3f);
        K2x3=h*(x1f.*x2f - beta.*x3f - dab3 + (p3f.*(miu.*x3f.^2 + 1))./2 + miu.*x3f.*(kapa - 1./2));
        
        x1(:,k+1)=x01(:,k)+1/2*(K1x1+K2x1);
        x2(:,k+1)=x02(:,k)+1/2*(K1x2+K2x2);
        x3(:,k+1)=x03(:,k)+1/2*(K1x3+K2x3);
    end
    
    K=length(x1(:,1));
    I=[];
    for k=1:K
        if (max(abs(x1(k,:)-x01(k,:))+abs(x2(k,:)-x02(k,:))+abs(x3(k,:)-x03(k,:)))<delta)&&(i>1000)
            xT1(pos(k))=x1(k,end);
            xT2(pos(k))=x2(k,end);
            xT3(pos(k))=x3(k,end);
            I=[I k];
        end
        if (max(abs(x1(k,:)-x01(k,:))+abs(x2(k,:)-x02(k,:))+abs(x3(k,:)-x03(k,:)))>10000)
            I0=[I0 pos(k)];
            I=[I k];
        end
    end
    if ~isempty(I)
        pos(I)=[];
        x1(I,:)=[];
        x2(I,:)=[];
        x3(I,:)=[];
        p1(I,:)=[];
        p2(I,:)=[];
        p3(I,:)=[];
    end
    
    if isempty(x1)
        break;
    end
    x01=x1;
    x02=x2;
    x03=x3;
    p01=p1;
    p02=p2;
    p03=p3;
end
xT1(I0)=[];
xT2(I0)=[];
xT3(I0)=[];
lambda1(I0)=[];
lambda2(I0)=[];
lambda3(I0)=[];

figure;
plot3(xT1,xT2,xT3,'*');

figure;
plot3(lambda1,lambda2,lambda3,'.');

figure;
plot(xT1,xT2,'*');
figure;
plot(xT1,xT3,'*');
figure;
plot(xT2,xT3,'*');
% 
% figure;
% plot(lambda1,lambda2,'.');
% figure;
% plot(lambda1,lambda3,'.');

% 
% m=m-length(pos);
% xT1(pos)=[];
% xT2(pos)=[];
% lambda1(pos)=[];
% lambda2(pos)=[];

%%% Neural Network
m=length(xT1);
eta=0.01;
A0=[xT1';xT2';xT3'];
Lambda=[lambda1';lambda2';lambda3'];
n1=20;
n2=20;
n3=20;
nf=3;
DELTAw=0.1;
DELTAb=0.01;
W1=DELTAw*randn(n1,3);
b1=DELTAb*randn(n1,1);
W2=DELTAw*randn(n2,n1);
b2=DELTAb*randn(n2,1);
W3=DELTAw*randn(n3,n2);
b3=DELTAb*randn(n3,1);
Wf=DELTAw*randn(nf,n3);
bf=DELTAb*randn(nf,1);
M=1/m*ones(m,1);
I=ones(1,m);
Jcost=[];

for i=1:500000
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
%     Af=20*tanh(Zf);
    Af=Zf;
%     I0=Af<-20;
%     Af(I0)=0;
%     I0=A2<0;
%     A2(I0)=0;
    
    Jcost=[Jcost sum((Af(1,:)-Lambda(1,:)).^2+(Af(2,:)-Lambda(2,:)).^2+(Af(3,:)-Lambda(3,:)).^2)/2/m];
    
    Ldaf=Af-Lambda;
%     Ldzf=Ldaf.*(1-tanh(Zf).^2)*20;
%     Ldzf=Ldaf.*(Zf>=-20);
    Ldzf=Ldaf;
%     Ldz2=Lda2.*(Z2>=0);
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
    
    W1=W1-eta*Jdw1;
    b1=b1-eta*Jdb1;
    W2=W2-eta*Jdw2;
    b2=b2-eta*Jdb2;
    W3=W3-eta*Jdw3;
    b3=b3-eta*Jdb3;
    Wf=Wf-eta*Jdwf;
    bf=bf-eta*Jdbf;
    
    if max(max(abs(Jdw1)))+max(abs(Jdb1))+max(max(abs(Jdw2)))+max(abs(Jdb2))+max(max(abs(Jdw3)))+max(abs(Jdb3))+max(max(abs(Jdwf)))+max(abs(Jdbf))<0.0001
%     if max(max(abs(Jdw1)))+max(abs(Jdb1))+max(max(abs(Jdw3)))+max(abs(Jdb3))+max(max(abs(Jdwf)))+max(abs(Jdbf))<0.001
        break;
    end
end

figure;
plot(Jcost);

%%% Test
xend=[1;1;1];
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
% af=20*tanh(zf);
af=zf;
% I0=af<-20;
% af(I0)=0;
% I0=a2<0;
% a2(I0)=0;

% Npath=1e3;
x1=xnode(1)+0.1*T;
x2=xnode(2)+0.1*T;
x3=xnode(3)+0.1*T;
p1=zeros(1,Npath+1);
p2=zeros(1,Npath+1);
p3=zeros(1,Npath+1);
p1(end)=af(1);
p2(end)=af(2);
p3(end)=af(3);
x01=x1;
x02=x2;
x03=x3;
p01=p1;
p02=p2;
p03=p3;
delta=0.0001;

for i=1:100000
    for k=1:Npath
        x1f=x01(:,Npath+2-k);
        x2f=x02(:,Npath+2-k);
        x3f=x03(:,Npath+2-k);
        p1f=p01(:,Npath+2-k);
        p2f=p02(:,Npath+2-k);
        p3f=p03(:,Npath+2-k);
        K1p1=h*(p3f.*x2f - p1f.*sigma + p2f.*(rou - x3f) + (miu.*x2f.*x3f)./(miu.*x3f.^2 + 1));
        K1p2=h*(p1f.*sigma - p2f + p3f.*x1f + (miu.*x1f.*x3f)./(miu.*x3f.^2 + 1));
        K1p3=h*((miu.*p3f.^2.*x3f)./2 - p2f.*x1f - p3f.*(beta - miu.*(kapa - 1./2)) +...
            (miu.*(x1f.*x2f - beta.*x3f + miu.*x3f.*(kapa - 1./2)))./(miu.*x3f.^2 + 1) -...
            (miu.*x3f.*(beta - miu.*(kapa - 1./2)))./(miu.*x3f.^2 + 1) -...
            (2.*miu.^2.*x3f.^2.*(x1f.*x2f - beta.*x3f + miu.*x3f.*(kapa - 1./2)))./(miu.*x3f.^2 + 1).^2 -...
            (2.*dab3.^2.*miu.*x3f)./(miu.*x3f.^2 + 1).^2);
        
        x1f=x01(:,Npath+1-k);
        x2f=x02(:,Npath+1-k);
        x3f=x03(:,Npath+1-k);
        p1f=p01(:,Npath+2-k)+K1p1;
        p2f=p02(:,Npath+2-k)+K1p2;
        p3f=p03(:,Npath+2-k)+K1p3;
        K2p1=h*(p3f.*x2f - p1f.*sigma + p2f.*(rou - x3f) + (miu.*x2f.*x3f)./(miu.*x3f.^2 + 1));
        K2p2=h*(p1f.*sigma - p2f + p3f.*x1f + (miu.*x1f.*x3f)./(miu.*x3f.^2 + 1));
        K2p3=h*((miu.*p3f.^2.*x3f)./2 - p2f.*x1f - p3f.*(beta - miu.*(kapa - 1./2)) +...
            (miu.*(x1f.*x2f - beta.*x3f + miu.*x3f.*(kapa - 1./2)))./(miu.*x3f.^2 + 1) -...
            (miu.*x3f.*(beta - miu.*(kapa - 1./2)))./(miu.*x3f.^2 + 1) -...
            (2.*miu.^2.*x3f.^2.*(x1f.*x2f - beta.*x3f + miu.*x3f.*(kapa - 1./2)))./(miu.*x3f.^2 + 1).^2 -...
            (2.*dab3.^2.*miu.*x3f)./(miu.*x3f.^2 + 1).^2);
        
        p1(:,Npath+1-k)=p01(:,Npath+2-k)+1/2*(K1p1+K2p1);
        p2(:,Npath+1-k)=p02(:,Npath+2-k)+1/2*(K1p2+K2p2);
        p3(:,Npath+1-k)=p03(:,Npath+2-k)+1/2*(K1p3+K2p3);
    end
    
    for k=1:Npath
        x1f=x01(:,k);
        x2f=x02(:,k);
        x3f=x03(:,k);
        p1f=p01(:,k);
        p2f=p02(:,k);
        p3f=p03(:,k);
        K1x1=h*(p1f./2 - dab1 - sigma.*(x1f - x2f));
        K1x2=h*(p2f./2 - dab2 - x2f + rou.*x1f - x1f.*x3f);
        K1x3=h*(x1f.*x2f - beta.*x3f - dab3 + (p3f.*(miu.*x3f.^2 + 1))./2 + miu.*x3f.*(kapa - 1./2));
        
        x1f=x01(:,k)+K1x1;
        x2f=x02(:,k)+K1x2;
        x3f=x03(:,k)+K1x3;
        p1f=p01(:,k+1);
        p2f=p02(:,k+1);
        p3f=p03(:,k+1);
        K2x1=h*(p1f./2 - dab1 - sigma.*(x1f - x2f));
        K2x2=h*(p2f./2 - dab2 - x2f + rou.*x1f - x1f.*x3f);
        K2x3=h*(x1f.*x2f - beta.*x3f - dab3 + (p3f.*(miu.*x3f.^2 + 1))./2 + miu.*x3f.*(kapa - 1./2));
        
        x1(:,k+1)=x01(:,k)+1/2*(K1x1+K2x1);
        x2(:,k+1)=x02(:,k)+1/2*(K1x2+K2x2);
        x3(:,k+1)=x03(:,k)+1/2*(K1x3+K2x3);
    end
    
    if (max(abs(x1-x01)+abs(x2-x02)+abs(x3-x03))<delta)&&(i>1000)&&(max(abs(x1-x01)+abs(x2-x02)+abs(x3-x03))>10000)
        break;
    end
    x01=x1;
    x02=x2;
    x03=x3;
    p01=p1;
    p02=p2;
    p03=p3;
end

figure;
plot(T,x1);
hold on
plot(T,x2);
hold on
plot(T,x3);
hold off

figure;
plot(T,p1);
hold on
plot(T,p2);
hold on
plot(T,p3);
hold off

figure;
plot3(x1,x2,x3);

%%% Validation
X=zeros(6,Npath+1);
X(:,1)=[x1(1);x2(1);x3(1);p1(1);p2(2);p3(2)];
t0=0;
for i=1:Npath
    X(:,i+1)=rk4(t0,h,X(:,i));
end

figure;
plot(T,X(1,:));
hold on
plot(T,X(2,:));
hold on
plot(T,X(3,:));
hold off

figure;
plot(T,X(4,:));
hold on
plot(T,X(5,:));
hold on
plot(T,X(6,:));
hold off

figure;
plot3(X(1,:),X(2,:),X(3,:));

