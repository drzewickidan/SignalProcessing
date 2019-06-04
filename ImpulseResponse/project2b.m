close all, clear all, clc
%% B1
R = 1000;
C = 0.000001;

% matlab
t = 0:0.000001:0.020;
x1 = (1/(R*C))*exp(-t/(R*C)).*(t>=0);
x2 = heaviside(t) + -1 * heaviside(t-0.01);
y = (conv(x1,x2)) * 0.000001;
y = 1 - y(1:20001);
y(10000:20001) = y(10000:20001) - 1;
plot(t, y, 'Linewidth', 2)
hold on

% simulated
data = xlsread('multisim_data_r.xlsx');
t = data(:,1);
y = data(:,2);
plot(t, y, '--r', 'Linewidth' , 2)

% measured
t = 0:0.00002:0.020;
data = xlsread('mydaq_data_r.xlsx');
plot(t, data, ':g', 'Linewidth', 2)


title('Unit Step Reponse of 1k Ohm Resistor')
xlabel('time (s)')
ylabel('voltage (V)')
legend('matlab', 'simulated', 'measured')
grid on

%% B2
R = 1000;
C = 0.000001;
A = 1;

figure
% ideal capacitor voltage
T = 0.010;
T_end = 0.020;
[V_c, ~] = rc_voltages(R, C, A, T, T_end); 
plot(t, V_c(t), 'Linewidth', 2) 
hold on

% discrete convolution
t = 0:0.000001:0.020;
x1 = (1/(R*C))*exp(-t/(R*C)).*(t>=0);
x2 = heaviside(t) + -1 * heaviside(t-0.01);
y = (conv(x1,x2)) * 0.000001;
plot(t, y(1:20001), '--y', 'Linewidth', 2);
hold on

% simulated
data = xlsread('multisim_data_c.xlsx');
t = data(:,1);
y = data(:,2);
plot(t, y, '--r', 'Linewidth' , 2)

% measured
t = 0:0.00002:0.020;
data = xlsread('mydaq_data_c.xlsx');
plot(t, data, ':g', 'Linewidth', 2)


title('Unit Step Reponse of 1 uF Capacitor')
xlabel('time (s)')
ylabel('voltage (V)')
legend('ideal', 'matlab' ,'simulated', 'measured')
grid on

%% B3
figure
% measured
t = 0:0.00002:0.020;
data = xlsread('mydaq_data_c.xlsx');
plot(t,data, 'Linewidth' , 2)
hold on

% fitted
t_o = 0.01;
V_o = 1;
V_c = @(tau, t) (t < t_o & t >= 0) .* (V_o .* (1 - exp(-t/(tau)))) ...
    + (t >= t_o) .* (V_o .* (1 - exp(-t_o/(tau))) .* exp(-(t - t_o)/(tau)));
yi=data(end);
idx=max(find(abs(y - yi) >= 0.37 * yi));
tau_est = data(idx);
tau_est = lsqcurvefit(V_c, tau_est, t', data);
plot(t, V_c(tau_est, t), '--', 'Linewidth', 2)
hold on

title('Unit Step Reponse of 1 uF Capacitor')
xlabel('time (s)')
ylabel('voltage (V)')
ylim([0 1.2])
legend('measured', 'fitted')
grid on

figure
% ideal impulse response
syms T t
f = @(T, t) (V_o .* (1 - exp(-t/(T))));
H_c = matlabFunction(diff(f(T, t)));

tau = R*C;
t = 0:0.000001:0.020;
plot(t, H_c(tau, t), 'Linewidth', 2)
hold on

% estimated impulse response
syms T t
f = @(T, t) (V_o .* (1 - exp(-t/(T))));
H_c = matlabFunction(diff(f(T, t)));

t = 0:0.000001:0.020;
plot(t, H_c(tau_est, t), '--' ,'Linewidth', 2)

title('Impulse Reponse of 1 uF Capacitor')
xlabel('time (s)')
ylabel('voltage rate of change (V/s)')
legend('ideal', 'estimated')
grid on