close all;
%% Problem 1-1
figure
f = @(t) (t-1)*(t>=1)-(t-2)*(t>=2)-(t-3)*(t>=3)+(t-4)*(t>=4);
fplot(f, [0 5])
ylim([0 1.1])
set(get(gca,'Children'), 'Linewidth',2)
grid on
title('Problem 1')
%% Problem 1-2
figure
f = @(t) (t-1)*(t>=1)-2*(t>=2)-1*(t-3)*(t>=3);
fplot(f, [-1 5])
ylim([-1.1 1.1])
set(get(gca,'Children'), 'Linewidth',2)
grid on
title('Problem 2')
%% Problem 1-3
figure
f = @(t) 5*(t)*(t>=0)-5*(t-1)*(t>=1)-5*(t>=1)+5*(t-1)*(t>=1)-5*(t>=2)-5*(t-2)*(t>=2);
fplot(f, [-5 5])
ylim([-0.1 5.1])
set(get(gca,'Children'), 'Linewidth',2)
grid on
title('Problem 3')
%% Problem 1-4
figure
f = @(t) (t)*(t>=0)+(t>=1)-(t>=2)-2*(t-2)*(t>=2)-(t>=3)+(t-3)*(t>=3);
fplot(f, [-1 5], 'linewidth', 2)
ylim([-0.1 5.1])
set(get(gca,'Children'), 'Linewidth',2)
grid on
title('Problem 4')
%% Problem 1-5
figure
f = @(t) 2*(t>=-1)-2*(t>=3)-0.5*(t)*(t>=0)+1*(t-1)*(t>=1)-0.5*(t-2)*(t>=2);
fplot(f, [-5 5])
ylim([-0.1 2.5])
set(get(gca,'Children'), 'Linewidth',2)
grid on
title('Problem 5')
%% Problem 2A
% Linearity
a=2;
b=3;
u=@(t) t>=0;
x1=@(t) u(t)-u(t-2);
x2=@(t) u(t)-u(t-3);
y1=@(t) a*(x1(t)+x1(t-1));
y2=@(t) b*(x2(t)+x2(t-1));
y_combined=@(t) y1(t)+y2(t);
x_combined=@(t) a*x1(t)+b*x2(t);
y=@(t) x_combined(t)+x_combined(t-1);
figure
subplot(2,1,1)
fplot(y_combined,[-5,5])
ylim([0 11])
set(get(gca,'Children'),'Linewidth',2)
xlabel('t')
title('A: Sum of Individual Outputs')
subplot(2,1,2)
fplot(y,[-5,5])
ylim([0 11])
set(get(gca,'Children'),'Linewidth',2)
xlabel('t')
title('A: Output of Summed Inputs')

% TV
u=@(t) t>0;
x=@(t) u(t)-u(t-5); 
y=@(t) x(t)+x(t-1); % INPUT
figure
subplot(2,1,1)
fplot(x,[-10 10])
set(get(gca,'Children'), 'Linewidth',2)
xlabel('t')
title('x(t)')
ylim([0 1.2])
subplot(2,1,2)
fplot(y,[-10 10], 'linewidth',2) 
set(get(gca,'Children'),'Linewidth',2)
xlabel('t')
title('A: y(t)')
y_delayed=@(t) y(t-3); % Delay y by 3
x_delayed=@(t) x(t-3); % Delay input by 3 and then put through system
y_delayed1=@(t) 3*x_delayed(t)+2*cos(pi*t/3);
figure
subplot(2,1,1)
fplot(y_delayed,[-10 10])
set(get(gca,'Children'), 'Linewidth',2)
xlabel('t')
title('A: y(t) delayed by 3')
subplot(2,1,2)
fplot(y_delayed1,[-10 10])
set(get(gca,'Children'), 'Linewidth',2)
xlabel('t')
title('A: x(t) delayed by 3 and input to system')


%% Problem 2B
a=2;
b=3;
u=@(t) t>=0;
x1=@(t) u(t)-u(t-2);
x2=@(t) u(t)-u(t-3);
y1=@(t) a*(3*x1(t)+2*cos(pi*t/3));
y2=@(t) b*(3*x2(t)+2*cos(pi*t/3));
y_combined=@(t) y1(t)+y2(t);
x_combined=@(t) a*x1(t)+b*x2(t);
y=@(t) 3*x_combined(t)+2*cos(pi*t/3);
figure
subplot(2,1,1)
fplot(y_combined,[-5,5])
ylim([-11 26])
set(get(gca,'Children'),'Linewidth',2)
xlabel('t')
title('B: Sum of Individual Outputs')
subplot(2,1,2)
fplot(y,[-5,5])
ylim([-11 21])
set(get(gca,'Children'),'Linewidth',2)
xlabel('t')
title('B: Output of Summed Inputs')

% TV
u=@(t) t>0;
x=@(t) u(t)-u(t-5); 
y=@(t) 3*x1(t)+2*cos(pi*t/3); % INPUT
figure
subplot(2,1,1)
fplot(x,[-10 10])
set(get(gca,'Children'), 'Linewidth',2)
xlabel('t')
title('x(t)')
ylim([0 1.2])
subplot(2,1,2)
fplot(y,[-10 10], 'linewidth',2) 
set(get(gca,'Children'),'Linewidth',2)
xlabel('t')
title('B: y(t) A')
y_delayed=@(t) y(t-3); % Delay y by 3
x_delayed=@(t) x(t-3); % Delay input by 3 and then put through system
y_delayed1=@(t) x_delayed(t)+x_delayed(t-1);
figure
subplot(2,1,1)
fplot(y_delayed,[-10 10])
set(get(gca,'Children'), 'Linewidth',2)
xlabel('t')
title('B: y(t) delayed by 3')
subplot(2,1,2)
fplot(y_delayed1,[-10 10])
set(get(gca,'Children'), 'Linewidth',2)
xlabel('t')
title('B: x(t) delayed by 3 and input to system')
%% Problem 2C
a=2;
b=3;
u=@(t) t>=0;
x1=@(t) u(t)-u(t-2);
x2=@(t) u(t)-u(t-3);
y1=@(t) a*(x1(-t));
y2=@(t) b*(x2(-t));
y_combined=@(t) y1(t)+y2(t);
x_combined=@(t) a*x1(t)+b*x2(t);
y=@(t) x_combined(-t);
figure
subplot(2,1,1)
fplot(y_combined,[-5,5])
ylim([-1 10])
set(get(gca,'Children'),'Linewidth',2)
xlabel('t')
title('C: Sum of Individual Outputs')
subplot(2,1,2)
fplot(y,[-5,5],'linewidth',2)
ylim([-1 10])
set(get(gca,'Children'),'Linewidth',2)
xlabel('t')
title('C: Output of Summed Inputs')

% TV
u=@(t) t>0;
x=@(t) u(t)-u(t-5); 
y=@(t) x(-t); % INPUT
figure
subplot(2,1,1)
fplot(x,[-10 10])
set(get(gca,'Children'), 'Linewidth',2)
xlabel('t')
title('C: x(t)')
ylim([0 1.2])
subplot(2,1,2)
fplot(y,[-10 10], 'linewidth',2) 
set(get(gca,'Children'),'Linewidth',2)
xlabel('t')
title('C: y(t)')
y_delayed=@(t) y(t-3); % Delay y by 3
x_delayed=@(t) x(t-3); % Delay input by 3 and then put through system
y_delayed1=@(t) x_delayed(t)+x_delayed(t-1);
figure
subplot(2,1,1)
fplot(y_delayed,[-10 10])
set(get(gca,'Children'), 'Linewidth',2)
xlabel('t')
title('C: y(t) delayed by 3')
subplot(2,1,2)
fplot(y_delayed1,[-10 10])
set(get(gca,'Children'), 'Linewidth',2)
xlabel('t')
title('C: x(t) delayed by 3 and input to system')
%% Problem 2D
% Linearity
n=0:5;
x1=0.8.^n;
x2=cos(n);
a1=2;
a2=3;
z=(a1*x1+a2*x2); 
y1=(n.^2).*z; % y[n] = n^2 * x[n]
z1=(n.^2).*x1;
z2=(n.^2).*x2;
y2=a1*z1+a2*z2;
figure
subplot(2,1,1)
stem(n,y1,'filled')
set(get(gca,'Children'), 'Linewidth',1.1)
xlim([-0.1 5.1])
title('D: Transform of Sum')
subplot(2,1,2)
stem(n,y2,'filled')
set(get(gca,'Children'), 'Linewidth',1.1)
xlim([-0.1 5.1])
title('D: Sum of Transforms')

% TV
n=0:5;                     
x=n;                       
y=x.*n.^2;                 % Unshifted system response
figure                      
subplot(2,1,1)             
stem(n+2,y,'filled')       % Delay output by 2 units
set(get(gca,'Children'), 'Linewidth',1.1)
xlim([1.9 7.1])            
title('D: Delayed Output')

n=n+2;
x2=n; 
y2=x2.*n.^2;
subplot(2,1,2)
stem(n,y2,'filled')
set(get(gca,'Children'), 'Linewidth',1.1)
xlim([1.9 7.1])
title('D: Delayed Input')
%% Problem 2E
% Linearity
n=0:5;
x1=0.8.^n;
x2=cos(n);
a1=2;
a2=3;
z=(a1*x1+a2*x2); 
y1=z.^2; % y[n] = x^2 * x[n]
z1=x1.^2;
z2=x2.^2;
y2=a1*z1+a2*z2;
figure
subplot(2,1,1)
stem(n,y1,'filled')
set(get(gca,'Children'), 'Linewidth',1.1)
xlim([-0.1 5.1])
title('E: Transform of Sum')
subplot(2,1,2)
stem(n,y2,'filled')
set(get(gca,'Children'), 'Linewidth',1.1)
xlim([-0.1 5.1])
title('E: Sum of Transforms')

% TV
n=0:5;                     
x=n;                       
y=x.^2;     % Unshifted system response
figure                      
subplot(2,1,1)             
stem(n+2,y,'filled')    % Delay output by 2 units
set(get(gca,'Children'), 'Linewidth',1.1)
xlim([1.9 7.1])            
title('E: Delayed Output')

n=n+2;
x2=n; 
y2=x2.^2;
subplot(2,1,2)
stem(n,y2,'filled')
set(get(gca,'Children'), 'Linewidth',1.1)
xlim([1.9 7.1])
title('E: Delayed Input')