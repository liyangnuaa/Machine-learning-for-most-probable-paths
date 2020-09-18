function y=fun(~,x)

global c dab1 dab2 dab3 sigma beta rou

x1f=x(1);
x2f=x(2);
x3f=x(3);
p1f=x(4);
p2f=x(5);
p3f=x(6);

y=zeros(size(x));
y(1)=2*p1f*c^2 - 2*dab1 - 2*sigma*(x1f - x2f);
y(2)=2*p2f*c^2 - 2*dab2 - 2*x2f + 2*rou*x1f - 2*x1f*x3f;
y(3)=2*p3f*c^2 - 2*dab3 - 2*beta*x3f + 2*x1f*x2f;
y(4)=-(2*p3f*x2f - 2*p1f*sigma + p2f*(2*rou - 2*x3f));
y(5)=-(2*p1f*sigma - 2*p2f + 2*p3f*x1f);
y(6)=-(- 2*beta*p3f - 2*p2f*x1f);