% COMM. COMPUTER ASSIGNMENT PART 1
% ZACK GOYETCHE
 
clear all
close all
clc
 
% TRIANGULAR PULSE PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
startTime_g1 = -16; % start at -25 seconds
endTime_g1 = 16; % end at 25 seconds
fs_g1 = 4; % number of samples per second
timeStep_g1 = 1/fs_g1; % size of time step between samples
numSamples_g1 = (endTime_g1-startTime_g1)*fs_g1; % total # samples
t_g1 = startTime_g1:timeStep_g1:endTime_g1; % time vector
t_g1Acf = startTime_g1:timeStep_g1/2:endTime_g1; % time vector for ACF
f_LLg1 = -fs_g1/2;
f_ULg1 = fs_g1/2;
freqStep_g1 = fs_g1/(numSamples_g1);
f_g1 = f_LLg1:freqStep_g1:f_ULg1; % frequency vector
f_g1Ifft = f_LLg1:freqStep_g1/2:f_ULg1; % frequency vector
normalize_g1 = fs_g1;
 
% RECTANGULAR PULSE PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
startTime_g2 = -16; % start at -3 seconds
endTime_g2 = 16; % end at 3 seconds
fs_g2 = 4; % number of samples per second
timeStep_g2 = 1/fs_g2; % size of time step between samples
numSamples_g2 = (endTime_g2-startTime_g2)*fs_g2; % total # samples
t_g2 = startTime_g2:timeStep_g2:endTime_g2; % time vector
t_g2Acf = startTime_g2:timeStep_g2/2:endTime_g2; % time vector for ACF
f_LLg2 = -fs_g2/2;
f_ULg2 = fs_g2/2;
freqStep_g2 = fs_g2/(numSamples_g2);
f_g2 = f_LLg2:freqStep_g2:f_ULg2; % frequency vector
f_g2Ifft = f_LLg2:freqStep_g2/2:f_ULg2; % frequency vector
normalize_g2 = fs_g2;
 
% CREATE TRIANGULAR PULSE ARRAY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i = 1:numSamples_g1+1
    if t_g1(i) < -2
        g1(i) = 0;
    end
    if t_g1(i) >= -2 && t_g1(i) <= 0
        g1(i) = 1/2*t_g1(i)+1;
    end
    if t_g1(i) > 0 && t_g1(i) <= 2
        g1(i) = -1/2*t_g1(i)+1;
    end
    
    if t_g1(i) > 2 
        g1(i) = 0;
    end
end
 
% CREATE RECTANGULAR PULSE ARRAY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i = 1:numSamples_g2+1
    if t_g2(i) < -0.5
        g2(i) = 0;
    end
    if t_g2(i) >= -0.5 && t_g2(i) <= 0.5
        g2(i) = 1;
    end
    if t_g2(i) > 0.5 
        g2(i) = 0;
    end
end
 
% PLOT - TIME/FREQ FOR RECT/TRI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
% G1 TRIANGULAR FUNCTION############################################
% time domain
figure(1)
stem(t_g1,g1),grid,grid minor
axis([-3 3 0 1.2])
% text(-0.5,1.25,'Energy = 1.333 J')
% title('Time Domain g1(t)')
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(01)g1_time_discrete','-dpng')
 
% freq domain
figure(2)
G1 = fftshift(fft(g1));
G1_scaled = (G1)/normalize_g1;
stem(f_g1,abs(G1_scaled)),grid,grid minor
axis([-1.5 1.5 0 2.25])
% text(-0.25,2.5,'Energy = 1.333 J')
% title('Frequency Domain G1(f)')
xlabel('Freqency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(02)g1_freq_discrete','-dpng')
 
% G2 RECTANGULAR FUNCTION###########################################
% time domain
figure(3)
stem(t_g2,g2),grid,grid minor
axis([-1.5 1.5 0 1.3])
% text(-0.25,1.25,'Energy = 1.000 J')
% title('Time Domain g2(t)')
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(03)g2_time_discrete','-dpng')
 
% freq domain
figure(4)
G2 = fftshift(fft(g2));
G2_scaled = (G2)/normalize_g2;
stem(f_g2,abs(G2_scaled)),grid,grid minor
axis([-2 2 0 1.3])
% text(-0.5,1.25,'Energy = 1.000 J')
% title('Frequency Domain G2(f)')
xlabel('Freqency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(04)g2_freq_discrete','-dpng')
 
% AUTO CORRELATION DIRECT METHOD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% triangular pulse
acfg1 = xcorr(g1);
% rectangular pulse
acfg2 = xcorr(g2);
 
% FIND MAX VAL FROM A.C. FOR NORMALIZATION~~~~~~~~~~~~~~~~~~~~~~~~~~
% triangular pulse
max_acfg1 = 0;
for i = 1:numel(acfg1)
    if acfg1(i) > max_acfg1
        max_acfg1 = acfg1(i);
    end
end
% rectangular pulse
max_acfg2 = 0;
for i = 1:numel(acfg2)
    if acfg2(i) > max_acfg2
        max_acfg2 = acfg2(i);
    end
end
 
% NORMALIZE THE AMPLITUDE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
acfg1_corrected = acfg1/max_acfg1;
acfg2_corrected = acfg2/max_acfg2;
 
% AUTO CORRELATION VIA ESD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
acfg1_ifft = abs(ifftshift(ifft((G1).^2)))/max_acfg1;
acfg2_ifft = abs(ifftshift(ifft((G2).^2)))/max_acfg2;
 
% PLOT THE AUTOCORRELATION FROM D.M. AND ESD METHOD~~~~~~~~~~~~~~~~~
 
% stem for g1 (triangular pulse)
figure(5)
stem(t_g1Acf,acfg1_corrected),grid,grid minor
% title('g1 Direct Method')
axis([-3 3 0 1.1])
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(05)g1_AC_DM_discrete','-dpng')
 
% stem for g2 (rectangular pulse)
figure(6)
% stem(f_g2,abs((G2/normalize_g2).^2));
stem(t_g2Acf,acfg2_corrected),grid,grid minor
% title('g2 Direct Method')
axis([-1.5 1.5 0 1.1])
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(06)g2_AC_DM_discrete','-dpng')
 
% stem for g1 (triangular pulse)
figure(7)
stem(t_g1./2,acfg1_ifft),grid,grid minor
% title('F^-^1\{\Psi_G_1(f)\}')
axis([-3 3 0 1.1])
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(07)g1_AC_ESD_discrete','-dpng')
 
% stem for g2 (rectangular pulse)
figure(8)
stem(t_g2./2,acfg2_ifft),grid,grid minor
str4 = 'Energy = ';
% title('F^-^1\{\Psi_G_2(f)\}')
axis([-1.5 1.5 0 1.1])
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(08)g2_AC_ESD_discrete','-dpng')
 
% CALCULATE ENERGY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% calculation done in time domain
e_g1 = trapz(t_g1,g1.^2)
e_g2 = trapz(t_g2,g2.^2)
 
% calculation done in freq domain
e_G1 = trapz(f_g1,G1.*conj(G1))/fs_g1^2
e_G2 = trapz(f_g2,G2.*conj(G2))/fs_g2^2
 
% CREATE IDEAL LOW PASS FILTERS ARRAY WITH BW = 2.0 HZ~~~~~~~~~~~~~~~
 
BW = 2;
LPF_LB = -BW;
LPF_UB = BW;
LPF1_Amp = 2; %H1 LPF amplitude
LPF2_Amp = 1; %H2 LPF amplitude
for i = 1:numel(f_g1)
    if f_g1(i) <= LPF_LB
        H1(i) = 0;
    end
    if f_g1(i) > LPF_LB && f_g1(i) < LPF_UB
        H1(i) = LPF1_Amp;
    end
    if f_g1(i) >= LPF_UB
        H1(i) = 0;
    end
end
 
for i = 1:numel(f_g2)
    if f_g2(i) <= LPF_LB
        H2(i) = 0;
    end
    if f_g2(i) > LPF_LB && f_g2(i) < LPF_UB
        H2(i) = LPF2_Amp;
    end
    if f_g2(i) >= LPF_UB 
        H2(i) = 0;
    end
end
 
% BRING LPFs INTO TIME DOMAIN
h1 = abs(ifftshift(ifft(H1)))*fs_g1;
h2 = abs(ifftshift(ifft(H2)))*fs_g2;
 
% PLOT OF LOW PASS FILTERS IN FREQ AND TIME DOMAIN~~~~~~~~~~~~~~~~~~
figure(9)
stem(f_g1,H1),grid,grid minor
axis([-3 3 0 LPF1_Amp+0.5])
% title('H1(f) - LPF - Freq Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(09)LPF1_freq_discrete','-dpng')
 
figure(10)
stem(f_g2,H2),grid,grid minor
axis([-3 3 0 LPF2_Amp+0.5])
% title('H2(f) - LPF - Freq Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(10)LPF2_freq_discrete','-dpng')
 
% FREQ DOMAIN SIGNAL MULTIPLIED BY FREQ DOMAIN LPF~~~~~~~~~~~~~~~~~~~
Y1 = G1_scaled.*H1;
Y2 = G2_scaled.*H2;
 
y1 = ifft(Y1);
y2 = ifft(Y2);
 
figure(11)
stem(t_g1,abs(y1)/2*fs_g1),grid,grid minor
axis([-3 3 0 1.25])
% title('g1(t) w/ LPF - Time Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(11)y1_time_discrete','-dpng')
 
figure(12)
stem(t_g2,abs(y2)*fs_g2),grid,grid minor
axis([-1.5 1.5 0 1.25])
% title('g2(t) w/ LPF - Time Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(12)y2_time_discrete','-dpng')
 
figure(13)
stem(f_g1,abs(Y1)/2),grid,grid minor
axis([-2 2 0 2.25])
% title('G1(f) w/ LPF - Freq Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(13)Y1_freq_discrete','-dpng')
 
figure(14)
stem(f_g2,abs(Y2)),grid,grid minor
axis([-2 2 0 1.3])
% title('G2(f) w/ LPF - Freq Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(14)Y2_freq_discrete','-dpng')
 
% PLOT OF y (y1 MIXED WITH y2)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Y = Y1.*Y2;
y = ifftshift(ifft(Y));
 
% stem of y (time domain)
figure(15)
stem(t_g1,abs((y)/2*fs_g1)),grid,grid minor
axis([-3 3 0 1.25])
xlabel('Time (s)')
ylabel('Amplitude')
%title('y(t) - Time Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(15)y_time_discrete','-dpng')
 
% stem of Y (frequency domain)
figure(16)
stem(f_g1,abs(Y)/2),grid,grid minor
axis([-1.5 1.5 0 2.75])
xlabel('Frequency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\disc\(16)y_freq_discrete','-dpng')
% title('Y(f) - Freq Domain')
 
% END OF FILE
