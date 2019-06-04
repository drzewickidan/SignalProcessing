clear, clc, close all
multisim_data_square = csvread('Squarewave_Multisim.csv');
multisim_data_sawtooth = csvread('Sawtooth_Multisim.csv');
%% C
N = 50;
A = 1*sqrt(2);
T = 0.01;

figure
[~, f_square] = fourier_analysis_square(A,T,N);
fplot(f_square, [0 5*T])
grid on; hold on

multisim_data_square_input = multisim_data_square(:,8);
t = linspace(0, 5*T, length(multisim_data_square_input));
plot(t, multisim_data_square_input)

legend('Approximation Signal', 'Multisim Simulation', 'Location', 'SouthOutside')
title('Fourier Analysis and Synthesis of Input Signal')
xlabel('time (s)')
ylabel('voltage (V)')
%% D
R = 1*10^3;    % 1 kOhm
C = 1*10^-6;    % 1 uF

figure
[tf_c,tf_r,w_c,w_r] = rc_bode_plot(R,C);
opts = bodeoptions;
opts.Title.FontSize = 14;
opts.FreqUnits = 'Hz';
opts.YLabel.String{1} = 'Gain';
opts.Title.String = '1 uF Capacitor Bode Plot';
bodeplot(tf_c, opts)

figure
opts.Title.String = '1 kOhm Resistor Bode Plot';
bodeplot(tf_r, opts)

fprintf('capacitor frequencies evaluted:\n')
disp(w_c)
fprintf('resistor frequencies evaluted:\n')
disp(w_r)




















%% E Capacitor 
R = 1e3;    % 1 kOhm
C = 1e-6;    % 1 uF
N = 50;
A = sqrt(2);
T = 0.01;
w = 2*pi/T;

[a_n, ~] = fourier_analysis_square(A,T,N);
multisim_data_square_c = multisim_data_square(:,5);

d_k = [] * ones(length(N));
for k=1:length(a_n)
    h_c = 1/(1+1j*(w)*(k-1)*R*C);
    term = h_c*a_n(k);
    d_k(k) = term;
end

f = @(t) 0;
for n = 0:1:N-1
    f = @(t) f(t) + A * d_k(n+1) * exp(1j*n*t*w) ;
end

figure
t = linspace(0, 0.05, 10001);
plot(t, f(t)) % matlab
hold on; grid on
t = linspace(0, 0.05, length(multisim_data_square_c));
plot(t, multisim_data_square_c) % multisim















%% E Resistor
R = 1*10^3;    % 1 kOhm
C = 1*10^-6;    % 1 uF
N = 50;
A = sqrt(2);
T = 0.01;
w = 2*pi/T;

[a_n, ~] = fourier_analysis_square(A,T,N);
multisim_data_square_r = multisim_data_square(:,2);
mydaq_data_square_r = xlsread('SquareWave_Resistor_Circuit.xlsx');

d_k = [] * ones(length(N));
for k=1:length(a_n)
    h_r = (1j*(w)*(k-1)*R*C)/(1+1j*(w)*(k-1)*R*C);
    term = h_r*a_n(k);
    d_k = [d_k term];
end

f = @(t) 0;
for n = 0:N-1
    f = @(t) f(t) + sqrt(2) * d_k(n+1) * exp(1j*n*t*w);
end

figure
t = linspace(0, 0.05, 10001);
plot(t, f(t)) % matlab
hold on; grid on

t = linspace(0, 0.05, length(multisim_data_square_r));
plot(t, multisim_data_square_r); % multisim

t = linspace(0, 0.05, length(mydaq_data_square_r));
plot(t, mydaq_data_square_r); % mydaq




































%% F-C
clear, clc, close all
N = 50;
A = 2;
T = 0.01;


[a_n, f_sawtooth] = fourier_analysis_sawtooth(A,T,N);

figure
fplot(f_sawtooth, [0 5*T]) % matlab
grid on; hold on

t = linspace(0, 5*T, length(multisim_data_sawtooth_input));
plot(t, multisim_data_sawtooth_input)

legend('Approximation Signal', 'Multisim Simulation', 'Location', 'SouthOutside')
title('Fourier Analysis and Synthesis of Input Signal')
xlabel('time (s)')
ylabel('voltage (V)')















%% F-D
R = 1*10^3;    % 1 kOhm
C = 1*10^-6;    % 1 uF

[tf_c,tf_r,w_c,w_r] = rc_bode_plot(1000,0.000001);

figure
opts = bodeoptions;
opts.Title.FontSize = 14;
opts.FreqUnits = 'Hz';
opts.YLabel.String{1} = 'Gain';
opts.Title.String = '1 uF Capacitor Bode Plot';
bodeplot(tf_c, opts)

figure
opts.Title.String = '1 kOhm Resistor Bode Plot';
bodeplot(tf_r, opts)

fprintf('capacitor frequencies evaluted:\n')
disp(w_c)
fprintf('resistor frequencies evaluted:\n')
disp(w_r)






















%% F-E Capacitor 
R = 1e3;    % 1 kOhm
C = 1e-6;    % 1 uF
N = 50;
A = sqrt(2);
T = 0.01;
w = 2*pi/T;

[a_n, ~] = fourier_analysis_sawtooth(A,T,N);
multisim_data_sawtooth_c = multisim_data_sawtooth(:,5);
mydaq_data_sawtooth_c = xlsread('Saw_Cap.xlsx');

d_k = [] * ones(length(N));
for k=1:length(a_n)
    h_c = (1j*(w)*(k-1)*R*C)/(1+1j*(w)*(k-1)*R*C);
    term = h_c*a_n(k);
    d_k(k) = term;
end

f = @(t) 0;
for n = 0:1:N-1
    f = @(t) f(t) + A * d_k(n+1) * exp(1j*n*t*w) ;
end

figure
t = linspace(0, 0.05, 10001);
%plot(t, f(t)) % matlab
hold on; grid on

t = linspace(0, 0.05, length(multisim_data_sawtooth_c));
plot(t, multisim_data_sawtooth_c) % multisim

t = linspace(0, 0.05, length(mydaq_data_sawtooth_c));
plot(t, mydaq_data_sawtooth_c)
legened('matlab', 'multisim', 'mydaq')
















%% F-E Resistor
R = 1*10^3;    % 1 kOhm
C = 1*10^-6;    % 1 uF
N = 50;
A = sqrt(2);
T = 0.01;
w = 2*pi/T;

[a_n, ~] = fourier_analysis_sawtooth(A,T,N);
mydaq_data_sawtooth_c = xlsread('Saw_Res.xlsx');



d_k = [] * ones(length(N));
for k=1:length(a_n)
    h_r = (1j*(w)*(k-1)*R*C)/(1+1j*(w)*(k-1)*R*C);
    term = h_r*a_n(k);
    d_k = [d_k term];
end

f = @(t) 0;
for n = 0:N-1
    f = @(t) f(t) + sqrt(2) * d_k(n+1) * exp(1j*n*t*w);
end

figure
t = linspace(0, 0.05, 10001);
plot(t, f(t))
hold on
t = linspace(0, 0.05, 5038);
plot(t, 10001)
legend('matlab', 'multisim')
t = linspace(0, 0.05, length(mydaq_data_sawtooth_c));
plot(t, mydaq_data_sawtooth_c)
