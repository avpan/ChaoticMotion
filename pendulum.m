% MCE 372 Engineering Analysis Example Code 
% Solution of Nonlinear Pendulum Problem 

function pendulum 

% Run Two Cases with Different Initial Positions 
clc;clear all;clf 
% g=9.81;m=100;L=1;Io=50; 
% K=m*g*L/Io; 

 %initial conditions [theta_0, w_0]
xo=[pi/2,0]; %initial conditions for theta0 and w0
[t,x]=ode45(@DE2,[0:0.01:5],xo); 
%figure(1);
%plot(x(:,1),x(:,2),'k--','Linewidth',1.5) 
%xlabel('\theta (radians)'),ylabel('Angular Velocity, w (rad/s)'),grid on 
%axis([0,5,-1.5,2]) 

%%plot radians vs time of Linear non-forced case
% xo=[pi/2,0];
% [t,x]=ode45(@DE2,[0:0.01:5],xo); 
% plot(t,x(:,1),'k-','Linewidth',1.5) 
% xlabel('t (1/k)'),ylabel('\theta  (radians)'),grid on 

% Plot Linearized Solution of Previous Two Cases 
h = 1E-3;
x1 = [pi/2+h,0];
xo = [pi/2,0]; 
[t,x]=ode45(@DE3,[0:.01:50],xo); 
[t,y]=ode45(@DE3,[0:.01:50],x1); 

%plot the chaotic pendlum
%figure(1);
%plot(t(:,1),x(:,1));

% figure(1);
% hold on;
% plot(x(:,1),x(:,2),'k--','Linewidth',1.5);
% plot(y(:,1),y(:,2),'k--','Linewidth',1.5, 'Color','r');
% xlabel('\theta (radians)'),ylabel('Angular Velocity, w (rad/s)'),grid on 

%remove t=0 comp, remove first row
X(:,1) = removerows(x(:,1),1);
X(:,2) = removerows(x(:,2),1);
Y(:,1) = removerows(y(:,1),1);
Y(:,2) = removerows(y(:,2),1);
T = removerows(t(:,1),1);

%find the Lyopunov exponent: dZ = dZ0exp(lambda*t)
dZ0 = sqrt(sum((x1-xo).^2));
dtheta = [Y(:,1) - X(:,1)];
dw = [Y(:,2) - X(:,2)];
dZ = sqrt(dtheta.^2 + dw.^2);
Lambda = log(dZ./dZ0)./T(:,1);
LambdaL= mean(Lambda)

for i=1:size(T(:,1))
    deltaZ(i,1) = (x1(1,1)-xo(1,1))*exp(LambdaL*T(i,1));
    deltaZ(i,2) = (x1(1,2)-xo(1,2))*exp(LambdaL*T(i,1));
end

figure(1);
semilogy(T(:,1),(deltaZ(:,1)),'k--','Linewidth',3);
xlabel('t (\omega_0^{-1})'),ylabel('\deltaZ')

clear all;
%x(1) is the theta 
%x(2) is the ang freq, w.
function dxdt=DE2(t,x) 
g=9.81;m=100;L=1.0; 
K = g/L;
dxdt=[x(2);-K*sin(x(1))]; 

function dxdt=DE3(t,x) 
%g=9.81;m=100;L=1.0;A =10.0;b=60;k=.5; 
% B = b/(m*L.^2)/k;
% K = g/L/k^2;
% f = A/k;
K = 3/2; q = 0;f = 3/4;
dxdt=[x(2);-K^2*sin(x(1))-q*x(2)+f*cos(t(1))]; 

% function dxdt=DE3(t,x) 
% g=9.81;m=100;L=1.0;A =10.0;b=60;k=.5; 
% B = b/(m*L.^2);
% K = g/L;
% dxdt=[x(2);A/(m*L.^2)*cos(k*t(1))-K*sin(x(1))-B*x(2)]; 
