function y=fun(~,x)

global dab1 dab2 dab3 sigma beta rou kapa miu

x1f=x(1);
x2f=x(2);
x3f=x(3);
p1f=x(4);
p2f=x(5);
p3f=x(6);

y=zeros(size(x));
y(1)=p1f/2 - dab1 - sigma*(x1f - x2f);
y(2)=p2f/2 - dab2 - x2f + rou*x1f - x1f*x3f;
y(3)=x1f*x2f - beta*x3f - dab3 + (p3f*(miu*x3f^2 + 1))/2 + miu*x3f*(kapa - 1/2);
y(4)=-(p3f*x2f - p1f*sigma + p2f*(rou - x3f) + (miu*x2f*x3f)/(miu*x3f^2 + 1));
y(5)=-(p1f*sigma - p2f + p3f*x1f + (miu*x1f*x3f)/(miu*x3f^2 + 1));
y(6)=-((miu*p3f^2*x3f)/2 - p2f*x1f - p3f*(beta - miu*(kapa - 1/2)) +...
    (miu*(x1f*x2f - beta*x3f + miu*x3f*(kapa - 1/2)))/(miu*x3f^2 + 1) -...
    (miu*x3f*(beta - miu*(kapa - 1/2)))/(miu*x3f^2 + 1) -...
    (2*miu^2*x3f^2*(x1f*x2f - beta*x3f + miu*x3f*(kapa - 1/2)))/(miu*x3f^2 + 1)^2 - (2*dab3^2*miu*x3f)/(miu*x3f^2 + 1)^2);
