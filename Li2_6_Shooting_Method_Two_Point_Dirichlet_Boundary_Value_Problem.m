% 打靶法求解迪利克雷边值问题
function [t,y]=shooting(p,q,f,tspan,x0f,varargin)
ga=x0f(1); gb=x0f(2);
f1=@(t,x)[x(2); -q(t)*x(1)-p(t)*x(2)];
f2=@(t,x)f1(t,x)+[0; f(t)];
[t,y1]=ode45(f1,tspan,[1;0],varargin{:});
[t,y2]=ode45(f1,tspan,[0;1],varargin{:});
[t,yp]=ode45(f2,tspan,[0;0],varargin{:});
m=(gb-ga*y1(end,1)-yp(end,1))/y2(end,1);
[t,y]=ode45(f2,tspan,[ga;m],varargin{:});
end

p=@(t)(1);
q=@(t)(1);
f=@(t)(-1);
x0=[0,0];
[t,y]=shooting(p,q,f,[0,1],x0);
x1=0:1/56:1;
y1=cos(x1)+(1-cos(1))/sin(1)*sin(x1)-1;
hold on;
plot(x1,y1,'or')
plot(t,y(:,1),'-b')
