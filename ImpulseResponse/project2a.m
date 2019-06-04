close all; clear; clc
%% PART A1
R = [600, 1000, 1200];
A = 5;
C = 1e-6;
T = 0.010;
T_end = 0.020;
t = 0:0.000001:0.020;

type('rc_voltages.m')

for m = 1:length(R)
    figure
    [V_c, V_r] = rc_voltages(R(m), C, A, T, T_end);
    plot(t, V_c(t), 'Linewidth', 2)
    hold on
    plot(t, V_r(t), 'Linewidth', 2)
    hold on
    plot(t, V_c(t) + V_r(t), '--c', 'Linewidth', 2)        
    xlabel('time (s)')
    ylabel('voltage (V)')
    legend('V_c', 'V_r', 'impulse response')
    title([num2str(R(m)), ' ohms'])
    grid on
end
%% PART A2
R = 1000;
A = 5;
C = 1e-6;
T = 0.010;
T_end = 0.020;
t = 0:0.000001:0.060;
delay = [0.007 0.012 0.017];
for i = 1:length(delay)
    figure,
    [V_c, V_r] = rc_voltages(R, C, A, T, T_end);
    plot(t, V_c(t) + V_c(t-delay(i)) + V_c(t-2*delay(i)), 'Linewidth', 2);
    hold on
    plot(t, V_r(t) + V_r(t-delay(i)) + V_r(t-2*delay(i)), 'Linewidth', 2)
    hold on
    plot(t, V_c(t) + V_r(t) + V_c(t-delay(i)) + V_r(t-delay(i)) + V_c(t-2*delay(i)) + V_r(t-2*delay(i)), '--c', 'Linewidth', 2)
    xlabel('time (s)')
    ylabel('voltage (V)')
    legend('V_c', 'V_r', 'impulse response')
    title([num2str(delay(i)*1000), ' ms'])
    grid on
end