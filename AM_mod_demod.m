close all
clear
clc

%% Generate modulating (message) signal
m = 1;
Am = 1;                                % Amplitude of modulating signal
fm = 100;                              % Frequency of modulating signal
fs = 100*fm;                           % Sampling frequency
Tm = 1/fm;                             % Time period of modulating signal
t = 0:1/fs:5*Tm;                       % Total time for simulation
ym = Am*sin(2*pi*fm*t);
figure, subplot(311), plot(t,ym), hold on, title ('Modulating Signal')


%% Generate carrier signal
Ac = Am/m;                             % Amplitude of carrier signal
fc = 20*fm;                            % Frequency of carrier signal
yc = Ac*sin(2*pi*fc*t);
subplot(312), plot(t,yc), hold on, ylabel ('Amplitude (Volts)'), title ('Carrier Signal')


%% Modulate using AM modulation
AM = Ac*sin(2*pi*fc*t).*(1+m*sin(2*pi*fm*t));
subplot(313), plot(t,AM,'r'), hold off, xlabel ('Time (s)'), title ('AM Signal')

L = length(t);
Nf = 2^ceil(log2(L));
AM_fft = fftshift(fft(AM,Nf));         % Frequency Response of AM Signal with expected tones at -fc and +fc
f = (-Nf/2:1:Nf/2-1)*fs/Nf;
figure, plot(f,abs(AM_fft)), title('Freq Response of AM'), xlabel('f(Hz)'), ylabel('|AM(F)|');


%% Demodulate the AM signal using Non-Coherent Detection
% Non-Coherent Detection Step 1: Envelope Detection
Vc = zeros(1,length(t));
Vc(1) = 0;                             % Initial capacitor voltage
for i = 2:length(AM)
    if AM(i) > Vc(i-1)                 % Diode on (charging)
        Vc(i) = AM(i);
    else                               % Diode off (discharging)
        Vc(i) = Vc(i-1) - 0.023*Vc(i-1);
    end
end

% Non-Coherent Detection Step 2: Low Pass RC Filter
fN = fs/2;
h = fir1(10,fc*2/fs, 'low');            
ym_rec = filter(h,1,Vc);
ym_rec = ym_rec - mean(ym_rec);

figure, subplot(211), plot(t, Vc);
title('Envelope detector output of AM signal'); xlabel('Time (s)'); ylabel('Amplitude (Volts)');
subplot(212), plot(t, ym_rec);
title('Non-coherent demodulated AM signal'); xlabel('Time (s)'); ylabel('Amplitude (Volts)');

ym_rec_fft = fftshift(fft(ym_rec,Nf));              % Frequency Response of retrieved message signal
f = (-Nf/2:1:Nf/2-1)*fs/Nf;
figure, plot(f,abs(ym_rec_fft)), title('Freq Response of demodulated message signal y_m(t)'), xlabel('f(Hz)'), ylabel('|AM(F)|');


%% Demodulate the AM signal using Coherent Detection
% Coherent Detection Step 1: Synchronous Demodulation using carrier
Vc = 2*AM.*sin(2*pi*fc*t);

% Coherent Detection Step 2: Low Pass RC Filter
[b,a] = butter(2,fc*2/fs);
ym_rec = filtfilt(b,a,Vc); % filtering the demodulated signal
ym_rec=ym_rec-mean(ym_rec);

figure, plot(t, ym_rec,'LineWidth',2), hold on, plot(t,ym,'r');
title('Coherent demodulated AM signal'); xlabel('Time (s)'); ylabel('Amplitude (Volts)');

ym_rec_fft = fftshift(fft(ym_rec,Nf));              % Frequency Response of retrieved message signal
f = (-Nf/2:1:Nf/2-1)*fs/Nf;
figure, plot(f,abs(ym_rec_fft)), title('Freq Response of demodulated message signal y_m(t)'), xlabel('f(Hz)'), ylabel('|AM(F)|');


