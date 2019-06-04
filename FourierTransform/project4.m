close all; clc; clear
%% 1uF
R = 1000;
L = 22 * 10^-3;
C = 1 * 10^-6;
a = R/(2*L);

figure
subplot(3,1,1)
hold on; 
k=sqrt(a^2-(1/(L*C)));
v = @(t) (2*a/k)*exp(-a.*t).*sinh(k.*t);
t=linspace(0, 0.001,1000);
plot(t,v(t))
t=xlsread('1uF_Multisim.xlsx','A:A');
resistorVoltage=xlsread('1uF_Multisim.xlsx','B:B');
plot(t,resistorVoltage)
t=xlsread('1uF_MyDaq.xlsx','A:A');
resistorVoltage=xlsread('1uF_MyDaq.xlsx','B:B');
plot(t,resistorVoltage)
legend('Theoretical','Simulated','Measured')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage Over Resistor with 1uF Capacitor')

subplot(3,1,2)
hold on
v = @(t) (1/k)*exp(-a.*t).*(k*cosh(k.*t)-a*sinh(k.*t));
t=linspace(0, 0.001,1000);
plot(t,v(t))
t=xlsread('1uF_Multisim.xlsx','D:D');
inductorVoltage=xlsread('1uF_Multisim.xlsx','E:E');
plot(t,inductorVoltage)
t=xlsread('1uF_MyDaq.xlsx','D:D');
inductorVoltage=xlsread('1uF_MyDaq.xlsx','E:E');
plot(t,inductorVoltage)
xlim([0,0.001])
legend('Theoretical','Simulated','Measured')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage Over Inductor with 1uF Capacitor')

subplot(3,1,3)
hold on
v = @(t) (1/(L*C*k*(a^2-k^2))).*(exp(-a.*t).*(-a.*sinh(k.*t)-k.*cosh(k.*t))+k);
t=linspace(0, 0.001,1000);
plot(t,v(t))
t=xlsread('1uF_Multisim.xlsx','G:G');
capacitorVoltage=xlsread('1uF_Multisim.xlsx','H:H');
plot(t,capacitorVoltage)
t=xlsread('1uF_MyDaq.xlsx','G:G');
capacitorVoltage=xlsread('1uF_MyDaq.xlsx','H:H');
plot(t,capacitorVoltage)
legend('Theoretical','Simulated','Measured')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage Over Capacitor with 1uF Capacitor')

%% 88nF
C = 88*10^-9;
a = R/(2*L);

figure
subplot(3,1,1)
hold on
v = @(t) 2.*a.*t.*exp(-a.*t);
t=linspace(0, 0.001,1000);
plot(t,v(t))
t=xlsread('88nF_Multisim.xlsx','A:A');
resistorVoltage=xlsread('88nF_Multisim.xlsx','B:B');
plot(t,resistorVoltage)
resistorVoltage=xlsread('88nF_MyDaq.xlsx','B:B');
t=xlsread('88nF_MyDaq.xlsx','A:A');
plot(t,resistorVoltage)
xlim([0,0.001])
legend('Theoretical','Simulated','Measured')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage Over Resistor with 88nF Capacitor')

subplot(3,1,2)
hold on
v = @(t) exp(-a.*t).*(1-a.*t);
t=linspace(0, 0.001,1000);
plot(t,v(t))
t=xlsread('88nF_Multisim.xlsx','D:D');
inductorVoltage=xlsread('88nF_Multisim.xlsx','E:E');
plot(t,inductorVoltage)
t=xlsread('88nF_MyDaq.xlsx','D:D');
inductorVoltage=xlsread('88nF_MyDaq.xlsx','E:E');
plot(t,inductorVoltage)
legend('Theoretical','Simulated','Measured')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage Over Inductor with 88nF Capacitor')

subplot(3,1,3)
hold on
v = @(t) (-1/(L*C))*(exp(-a.*t).*(1/a^2+t./a)+1/a^2)+2;
t=linspace(0, 0.001,1000);
plot(t,v(t))
t=xlsread('88nF_Multisim.xlsx','G:G');
capacitorVoltage=xlsread('88nF_Multisim.xlsx','H:H');
plot(t,capacitorVoltage)
capacitorVoltage=xlsread('88nF_MyDaq.xlsx','H:H');
t=xlsread('88nF_MyDaq.xlsx','G:G');
plot(t,capacitorVoltage)
legend('Theoretical','Simulated','Measured')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage Over Capacitor with 88nF Capacitor')

%% 300pF 
figure
C=300*10^-12;
w0=sqrt(1/(L*C)-a^2);
subplot(3,1,1)
xlim([0,0.001])
hold on
v = @(t) (2*a/w0)*exp(-a.*t).*sin(w0.*t);
t=linspace(0, 0.001,1000);
plot(t,v(t))
t=xlsread('300pF_Multisim.xlsx','A:A');
resistorVoltage=xlsread('300pF_Multisim.xlsx','B:B');
plot(t,resistorVoltage)
resistorVoltage=xlsread('300pF_MyDaq.xlsx','B:B');
t=xlsread('300pF_MyDaq.xlsx','A:A');
plot(t,resistorVoltage)
legend('Theoretical','Simulated','Measured')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage Over Resistor (300pF Capacitor)')

subplot(3,1,2)
hold on
v = @(t) (1/w0)*exp(-a.*t).*(w0*cos(w0.*t)-a*sin(w0.*t));
t=linspace(0, 0.001,1000);
plot(t,v(t))
t=xlsread('300pF_Multisim.xlsx','D:D');
inductorVoltage=xlsread('300pF_Multisim.xlsx','E:E');
plot(t,inductorVoltage)
inductorVoltage=xlsread('300pF_MyDaq.xlsx','E:E');
t=xlsread('300pF_MyDaq.xlsx','D:D');
plot(t,inductorVoltage)
legend('Theoretical','Simulated','Measured')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage Over Inductor (300pF Capacitor)')

subplot(3,1,3)
hold on
v = @(t) (1/(L*C*w0*(a^2+w0^2)))*(exp(-a.*t).*(a.*sin(w0.*t)-w0*cos(w0.*t))+w0);
t=linspace(0, 0.001,1000);
plot(t,v(t))
t=xlsread('300pF_Multisim.xlsx','G:G');
capacitorVoltage=xlsread('300pF_Multisim.xlsx','H:H');
plot(t,capacitorVoltage)
capacitorVoltage=xlsread('300pF_MyDaq.xlsx','H:H');
t=xlsread('300pF_MyDaq.xlsx','G:G');
plot(t,capacitorVoltage)
xlim([0,0.001])
legend('Theoretical','Simulated','Measured')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltage Over Capacitor (300pF Capacitor)')

%%
clear
w = logspace(0,8,5000);
R = 1e3; 
L = 22e-3; 
C = [300e-12, 88e-9, 1e-6];

Hr = @(w, C) 1j*w*R*C/(1-w^2*L*C+1j*w*R*C);
Hl = @(w, C) (-w^2*L*C)/(1-w^2*L*C+1j*w*R*C);
Hc = @(w, C) 1/(1-w^2*L*C+1j*w*R*C);

Res = [];
Ind = [];
Cap = [];
for c_i=1:length(C) 

    for i=1:length(w)
        Res = [Res Hr(w(i), C(c_i))];
        Ind = [Ind Hl(w(i), C(c_i))];
        Cap = [Cap Hc(w(i), C(c_i))];
    end
end

% 300pf
figure
subplot(2,1,2)
semilogx(w,rad2deg(angle(Res(1:5000))))
hold on
semilogx(w,rad2deg(angle(Ind(1:5000))))
hold on
semilogx(w,rad2deg(angle(Cap(1:5000))))

title('300pF Capacitor Phase Plot')
legend('H_R(w)','H_L(w)','H_C(w)')
xlabel('Frequency (w)')
ylabel('Phase (deg)')

subplot(2,1,1)
semilogx(w,mag2db(abs(Res(1:5000))))
hold on
semilogx(w,mag2db(abs(Ind(1:5000))))
hold on
semilogx(w,mag2db(abs(Cap(1:5000))))

title('300pf Capacitor Magnitutde Plot')
legend('H_R(w)','H_L(w)','H_C(w)')
xlabel('Frequency (w)')
ylabel('Magnitude (dB)')

% 1uF
figure
subplot(2,1,2)
semilogx(w,rad2deg(angle(Res(10001:15000))))
hold on
semilogx(w,rad2deg(angle(Ind(10001:15000))))
hold on
semilogx(w,rad2deg(angle(Cap(10001:15000))))

title('1uF Capacitor Phase Plot')
legend('H_R(w)','H_L(w)','H_C(w)')
xlabel('Frequency (w)')
ylabel('Phase (deg)')

subplot(2,1,1)
semilogx(w,mag2db(abs(Res(5001:10000))))
hold on
semilogx(w,mag2db(abs(Ind(5001:10000))))
hold on
semilogx(w,mag2db(abs(Cap(5001:10000))))

title('1uF Capacitor Magnitutde Plot')
legend('H_R(w)','H_L(w)','H_C(w)')
xlabel('Frequency (w)')
ylabel('Magnitude (dB)')

% 88nF
figure
subplot(2,1,2)
semilogx(w,rad2deg(angle(Res(5001:10000))))
hold on
semilogx(w,rad2deg(angle(Ind(5001:10000))))
hold on
semilogx(w,rad2deg(angle(Cap(5001:10000))))

title('88nF Capacitor Phase Plot')
legend('H_R(w)','H_L(w)','H_C(w)')
xlabel('Frequency (w)')
ylabel('Phase (deg)')

subplot(2,1,1)
semilogx(w,mag2db(abs(Res(10001:15000))))
hold on
semilogx(w,mag2db(abs(Ind(10001:15000))))
hold on
semilogx(w,mag2db(abs(Cap(10001:15000))))


title('88nF Capacitor Magnitutde Plot')
legend('H_R(w)','H_L(w)','H_C(w)')
xlabel('Frequency (w)')
ylabel('Magnitude (dB)')