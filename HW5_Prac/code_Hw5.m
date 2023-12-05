% Amirreza Hatamipour
% 97101507
%% Question 1
%x_(n+1)=Ax_n(1-x_n)
A=[0.1 0.5 1 1.5 2 2.5 3 3.5 4];
n=1;
for a=A
    x=zeros(1,500);
    x(1,1)=0.2;
    for i=1:499
        x(1,i+1)=a*x(1,i)*(1-x(1,i));
    end
    subplot(3,3,n)
    plot(x)
    title(a)
    grid on
    n=n+1;
end
%% section b
clear;clc; close all;
A=0:0.02:4;
n=1;
final=zeros(100,length(A));
for a=A
    x=zeros(1,500);
    x(1,1)=0.2;
    for i=1:499
        x(1,i+1)=a*x(1,i)*(1-x(1,i));
    end
    final(1:100,n)=x(1,401:500)';
    n=n+1;
end
for i=1:100
plot(A,final(i,:),'.')
hold on
title('bifurcartion')
grid on
end
%% Question 2
sigma=10;
b=2.67;
r=20.5;
[x,y,z] = lorenz(r, sigma, b);

subplot(3,1,1)
plot(x)
title('lorenz & r<24.74')
ylabel('x')
grid on
subplot(3,1,2)
plot(y)
ylabel('y')
grid on
subplot(3,1,3)
plot(z)
ylabel('z')
grid on
figure
plot3(x,y,z);
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('lorenz & r<24.74')
%% 24.74<r<28
r=26.5;
[x,y,z] = lorenz(r, sigma, b);

subplot(3,1,1)
plot(x)
title('lorenz & 24.74<r<28')
ylabel('x')
grid on
subplot(3,1,2)
plot(y)
ylabel('y')
grid on
subplot(3,1,3)
plot(z)
ylabel('z')
grid on

figure
plot3(x,y,z);
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('lorenz & 24.74<r<28')
%% r>28
r=40;
[x,y,z] = lorenz(r, sigma, b);

subplot(3,1,1)
plot(x)
title('lorenz & r>28')
ylabel('x')
grid on
subplot(3,1,2)
plot(y)
ylabel('y')
grid on
subplot(3,1,3)
plot(z)
ylabel('z')
grid on

figure
plot3(x,y,z);
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('lorenz & r>28')
%% Question 3
clc;clear; close all;

freq=[250 256 250 256 256 250 250 256 256 250];
feature_0=zeros(30,4);
for i=1:10
   filename="0_sig"+num2str(i)+"_"+num2str(freq(i))+".mat";
   data=load('seizure\'+filename);
   data=data.rec_signal;
   data1=data(1,1:3*freq(i));
   data2=data(1,1+freq(i):4*freq(i));
   data3=data(1,1+2*freq(i):5*freq(i));
   
   % feature one
   feature_0(1+(i-1)*3,1)=lyapunovExponent(data1,3*freq(i));
   feature_0(2+(i-1)*3,1)=lyapunovExponent(data2,3*freq(i));
   feature_0(3+(i-1)*3,1)=lyapunovExponent(data3,3*freq(i));
   
   %feature two
   feature_0(1+(i-1)*3,2)=approximateEntropy(data1);
   feature_0(2+(i-1)*3,2)=approximateEntropy(data2);
   feature_0(3+(i-1)*3,2)=approximateEntropy(data3);
   
   %feature three
   feature_0(1+(i-1)*3,3)=entropy(data1);
   feature_0(2+(i-1)*3,3)=entropy(data2);
   feature_0(3+(i-1)*3,3)=entropy(data3);
   
   %feature four
   feature_0(1+(i-1)*3,4)=correlationDimension(data1);
   feature_0(2+(i-1)*3,4)=correlationDimension(data2);
   feature_0(3+(i-1)*3,4)=correlationDimension(data3);
end

freq=[250 256 250 256 256 256 256 250 250 250];
feature_1=zeros(30,4);
for i=1:10
   filename="1_sig"+num2str(i)+"_"+num2str(freq(i))+".mat";
   data=load('seizure\'+filename);
   data=data.rec_signal;
   data1=data(1,1:3*freq(i));
   data2=data(1,1+freq(i):4*freq(i));
   data3=data(1,1+2*freq(i):5*freq(i));
   
   % feature one
   feature_1(1+(i-1)*3,1)=lyapunovExponent(data1,3*freq(i));
   feature_1(2+(i-1)*3,1)=lyapunovExponent(data2,3*freq(i));
   feature_1(3+(i-1)*3,1)=lyapunovExponent(data3,3*freq(i));
   
   %feature two
   feature_1(1+(i-1)*3,2)=approximateEntropy(data1);
   feature_1(2+(i-1)*3,2)=approximateEntropy(data2);
   feature_1(3+(i-1)*3,2)=approximateEntropy(data3);
   
   %feature three
   feature_1(1+(i-1)*3,3)=entropy(data1);
   feature_1(2+(i-1)*3,3)=entropy(data2);
   feature_1(3+(i-1)*3,3)=entropy(data3);
   
   %feature four
   feature_1(1+(i-1)*3,4)=correlationDimension(data1);
   feature_1(2+(i-1)*3,4)=correlationDimension(data2);
   feature_1(3+(i-1)*3,4)=correlationDimension(data3);
end

mean_feature_0=mean(feature_0)
var_feature_0=var(feature_0)

mean_feature_1=mean(feature_1)
var_feature_1=var(feature_1)

%% calculate var and mean after normalize feature
normalize_feature=mapminmax([feature_0;feature_1]')';
mean_normalize_feature_0=mean(normalize_feature(1:30,:))
var_normalize_feature_0=var(normalize_feature(1:30,:))

mean_normalize_feature_1=mean(normalize_feature(31:end,:))
var_normalize_feature_1=var(normalize_feature(31:end,:))   
%% Question 4

syms n
box_counting_cantor_set=vpa(limit(-log(2^n)/log((1/3)^n), n, inf, 'left'))
box_counting_koch_set=vpa(limit(-log(4^n)/log((1/3)^n), n, inf, 'left'))

%% Function
function [x,y,z] = lorenz(rho, sigma, beta, initV, T, eps)
% LORENZ Function generates the lorenz attractor of the prescribed values
% of parameters rho, sigma, beta
%
%   [X,Y,Z] = LORENZ(RHO,SIGMA,BETA,INITV,T,EPS)
%       X, Y, Z - output vectors of the strange attactor trajectories
%       RHO     - Rayleigh number
%       SIGMA   - Prandtl number
%       BETA    - parameter
%       INITV   - initial point
%       T       - time interval
%       EPS     - ode solver precision
%
% Example.
%        [X Y Z] = lorenz(28, 10, 8/3);
%        plot3(X,Y,Z);
if nargin<3
  error('MATLAB:lorenz:NotEnoughInputs','Not enough input arguments.'); 
end
if nargin<4
  eps = 0.000001;
  T = [0 25];
  initV = [0 1 1.05];
end
options = odeset('RelTol',eps,'AbsTol',[eps eps eps/10]);
[T,X] = ode45(@(T,X) F(T, X, sigma, rho, beta), T, initV, options);
plot3(X(:,1),X(:,2),X(:,3));
axis equal;
grid;
title('Lorenz attractor');
xlabel('X'); ylabel('Y'); zlabel('Z');
x = X(:,1);
y = X(:,2);
z = X(:,3);
return
end
function dx = F(T, X, sigma, rho, beta)
% Evaluates the right hand side of the Lorenz system
% x' = sigma*(y-x)
% y' = x*(rho - z) - y
% z' = x*y - beta*z
% typical values: rho = 28; sigma = 10; beta = 8/3;
    dx = zeros(3,1);
    dx(1) = sigma*(X(2) - X(1));
    dx(2) = X(1)*(rho - X(3)) - X(2);
    dx(3) = X(1)*X(2) - beta*X(3);
    return
end

