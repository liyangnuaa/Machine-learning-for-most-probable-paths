function y=fun(~,x)

global c dab1 Ch Gamma theta S0

x1f=x(1);
p1f=x(2);

y=zeros(size(x));
y(1)=((S0.*(tanh(x1f./10 - 53./2)./5 + 1./2))./4 - Gamma.*theta.*x1f.^4)./Ch - dab1 + (c.^2.*p1f)./2;
y(2)=-(- ((S0.*tanh(x1f./10 - 53./2).*(tanh(x1f./10 - 53./2).^2./10 - 1./10))./100 -...
     12.*Gamma.*theta.*x1f.^2)./Ch - (p1f.*((S0.*(tanh(x1f./10 - 53./2).^2./50 - 1./50))./4 + 4.*Gamma.*theta.*x1f.^3))./Ch);
