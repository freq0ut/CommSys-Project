% COMM. COMPUTER ASSIGNMENT PART 2
% ZACK GOYETCHE
 
close all
clear all
clc
 
% TRIANGULAR PULSE PARAMETERS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
startTime_g1 = -1; % start at -3 seconds
endTime_g1 = 1; % end at 3 seconds
fs_g1 = 4096; % number of samples per second
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
 
% CREATE TRIANGULAR PULSE ARRAY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i = 1:numSamples_g1+1
    if t_g1(i) < -0.04
        g1(i) = 0;
    end
    
    if t_g1(i) >= -0.04 && t_g1(i) < 0.02
        g1(i) = (2/0.02)*t_g1(i)+4;
    end
    
    if t_g1(i) >= -0.02 && t_g1(i) < 0
        g1(i) = -(2/0.02)*t_g1(i);
    end
    
    if t_g1(i) >= 0 && t_g1(i) < 0.02
        g1(i) = -(3/0.02)*t_g1(i);
    end
    
    if t_g1(i) >= 0.02 && t_g1(i) <= 0.04
        g1(i) = (3/0.02)*t_g1(i)-6;
    end
    
    if t_g1(i) > 0.04 
        g1(i) = 0;
    end
end
 
% PLOT - TIME/FREQ FOR RECT/TRI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(1)
% G1 TRIANGULAR FUNCTION############################################
% time domain
plot(t_g1,g1)
axis([-0.1 0.1 -3.5 3.5])
% title('Time Domain g1(t) [message signal]')
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(01)mt_timeDomain','-dpng')
 
% freq domain
figure(2)
G1 = fftshift(fft(g1));
G1_scaled = (G1)/normalize_g1; %scaled to 1
 
 
max_G1 = 0;
for i = 1:numel(G1_scaled)
    if abs(G1_scaled(i)) > max_G1
        max_G1 = abs(G1_scaled(i));
    end
end
plot(f_g1,abs(G1_scaled)),grid,grid minor
axis([-150 150 0 0.1]) %normalized
% title('Frequency Domain G1(f) [message signal]')
xlabel('Freqency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(02)mt_freqDomain','-dpng')
 
% %##################################################################
 
% CALCULATE ENERGY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%calculation done in time domain
e_g1 = trapz(t_g1,g1.^2)
 
% calculation done in freq domain
e_G1 = trapz(f_g1,G1.*conj(G1))/fs_g1^2
 
% DEFINE CARRIER~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f = 300;
c_g1 = cos(2*pi*f*t_g1);
 
C_g1 = fftshift(fft(c_g1));
C_scaled = (C_g1)/normalize_g1;
 
figure(3)
subplot(2,1,1)
plot(t_g1,c_g1)
axis([-1/60 1/60 -1.5 1.5])
title('c(t) Carrier Signal - Time Domain')
xlabel('Time (s)')
ylabel('Amplitude')
 
subplot(2,1,2)
plot(f_g1,abs(C_scaled))
axis([-350 350 0 1.1])
title('C(f) Carrier Signal - Frequency Domain')
xlabel('Freqency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(03)ct_time&freqDomain','-dpng')
 
% MIX MESSAGE AND CARRIER - PLOT IN TIME AND FREQ~~~~~~~~~~~~~~~~~~~
x = c_g1.*g1;
carrier_term = 3.5.*ones(1,numSamples_g1+1);
x_wc = c_g1.*(carrier_term+g1);
X = fftshift(fft(x))/normalize_g1;
X_wc = fftshift(fft(x_wc))/normalize_g1;
 
figure(4)
plot(t_g1,x)
axis([-0.05 0.05 -3.5 3.5])
xlabel('Time (s)')
ylabel('Amplitude')
% title('x(t) DSB-SC - Time Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(04)dsb-sc_timeDomain','-dpng')
 
figure(5)
plot(f_g1,abs(X))
axis([-450 450 0 0.06])
xlabel('Freqency (Hz)')
ylabel('Amplitude')
% title('X(f) DSB-SC - Freq Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(05)dsb-sc_freqDomain','-dpng')
 
figure(6)
plot(t_g1,x_wc)
axis([-0.05 0.05 -7.5 7.5])
xlabel('Time (s)')
ylabel('Amplitude')
% title('x_w_c(t) DSB-C - Time Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(06)dsb-c_timeDomain','-dpng')
 
figure(7)
plot(f_g1,abs(X_wc))
axis([-450 450 0 0.06])
xlabel('Freqency (Hz)')
ylabel('Amplitude')
% title('X_w_c(f) DSB-C - Freq Domain')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(07)dsb-c_freqDomain','-dpng')
 
% CREATE BUTTERWORTH FILTERS... 5TH and 40TH ORDER
% roll off at -(nth Order)*20dB/dec
 
n1=5; %nth order for filter #1
n2=40; %nth order for filter #2
ampH1 = 1;
ampH2 = 1;
 
rollOff_freq = 300; %cutoff frequency in Hz
plusDecade = 10*rollOff_freq; %used for determining slope of roll off
 
% FIND EQUATIONS FOR ROLL OFF LINES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ################################ H1 ####################################
% POSITIVE SLOPE
% num_slope_left H1
numSlopeLeft_H1 = -(10^n1) - ampH1;
% den_slope_left H1
denSlopeLeft_H1 = -plusDecade+rollOff_freq;
% combined left H1
slopeLeft_H1 = numSlopeLeft_H1/denSlopeLeft_H1;
% y intercept left H1
bH1_left = ampH1 - slopeLeft_H1*(-rollOff_freq);
% x intercept left H1
xIntercept_leftH1 = -bH1_left/slopeLeft_H1;
 
% NEGATIVE SLOPE
% num_slope_right H1
numSlopeRight_H1 = -(10^n1) - ampH1;
% den_slope_right H1
denSlopeRight_H1 = plusDecade-rollOff_freq;
% combined right H1
slopeRight_H1 = numSlopeRight_H1/denSlopeRight_H1;
% y intercept right H1
bH1_right = ampH1 - slopeRight_H1*(rollOff_freq);
% x intercept right H1
xIntercept_rightH1 = -bH1_right/slopeRight_H1;
 
%################################ H2 ######################################
% POSITIVE SLOPE
% num_slope_left H2
numSlopeLeft_H2 = -(10^n2) - ampH2;
% den_slope_left H2
denSlopeLeft_H2 = -plusDecade+rollOff_freq;
% combined left H2
slopeLeft_H2 = numSlopeLeft_H2/denSlopeLeft_H2;
% y intercept left H2
bH2_left = ampH2 - slopeLeft_H2*(-rollOff_freq);
% x intercept left H2
xIntercept_leftH2 = -bH2_left/slopeLeft_H2;
 
% NEGATIVE SLOPE
% num_slope_right H2
numSlopeRight_H2 = -(10^n2) - ampH2;
% den_slope_right H2
denSlopeRight_H2 = plusDecade-rollOff_freq;
% combined right H2
slopeRight_H2 = numSlopeRight_H2/denSlopeRight_H2;
% y intercept right H2
bH2_right = ampH2 - slopeRight_H2*(rollOff_freq);
% x intercept right H2
xIntercept_rightH2 = -bH2_right/slopeRight_H2;
 
% line H2 left: slopeLeft_H2*f_g1(i)+bH2_left
% line H2 right: slopeRight_H2*f_g1(i)+bH2_right
 
%DEFINE H1 AS A FILTER OF ORDER 5
for i = 1:numSamples_g1+1
    if f_g1(i) < xIntercept_leftH1
        H1(i) = 0;
    end
    
    if f_g1(i) >= xIntercept_leftH1 && f_g1(i) < -rollOff_freq
        H1(i) = slopeLeft_H1*f_g1(i)+bH1_left;
    end
    
    if f_g1(i) >= -rollOff_freq && f_g1(i) < rollOff_freq
        H1(i) = ampH1;
    end
    
    if f_g1(i) >= rollOff_freq && f_g1(i) < xIntercept_rightH1
        H1(i) = slopeRight_H1*f_g1(i)+bH1_right;
    end
    
    if f_g1(i) >= xIntercept_rightH1
        H1(i) = 0;
    end
end
 
for i = 1:numSamples_g1+1
    if f_g1(i) < xIntercept_leftH2
        H2(i) = 0;
    end
    
    if f_g1(i) >= xIntercept_leftH2 && f_g1(i) < -rollOff_freq
        H2(i) = slopeLeft_H2*f_g1(i)+bH2_left;
    end
    
    if f_g1(i) >= -rollOff_freq && f_g1(i) < rollOff_freq
        H2(i) = ampH2;
    end
    
    if f_g1(i) >= rollOff_freq && f_g1(i) < xIntercept_rightH2
        H2(i) = slopeRight_H2*f_g1(i)+bH2_right;
    end
    
    if f_g1(i) >= xIntercept_rightH2
        H2(i) = 0;
    end
end
 
figure(8)
plot(f_g1,H1)
axis([-600 600 0 1.5])
% title('Frequency Domain H1(f) [5^t^h Order]')
xlabel('Freqency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(08)5th_orderLPF','-dpng')
 
figure(9)
plot(f_g1,H2)
axis([-600 600 0 1.5])
% title('Frequency Domain H2(f) [40^t^h Order]')
xlabel('Freqency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(09)40th_orderLPF','-dpng')
 
h1 = ifft(H1);
h2 = ifft(H2);
 
demod_x = (x.*(2*c_g1));
demod_xplot = 4+x.*(2*c_g1);
demod_X = fftshift(fft(demod_x));
demod_Xplot = fftshift(fft(demod_xplot));
 
demod_x_wc = (x_wc.*(2*c_g1));
demod_X_wc = fftshift(fft(demod_x_wc));
 
Y_H1 = (demod_X).*H1;
Y_H2 = (demod_X).*H2;
 
Y_H1_plot = (demod_Xplot).*H1;
Y_H2_plot = (demod_Xplot).*H2;
 
Ywc_H1 = (demod_X_wc).*H1;
Ywc_H2 = (demod_X_wc).*H2;
 
% DEMOD DSB-SC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(10)
plot(t_g1,(demod_x))
axis([-0.1 0.1 -7 7])
% title('Demodulated Signal z(t)')
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(10)demodulated_dsb-sc_time','-dpng')
 
figure(11)
plot(f_g1,abs(demod_X)/normalize_g1)
axis([-700 700 0 0.1])
% title('Demodulated Signal Z(f)')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(11)demodulated_dsb-sc_freq','-dpng')
 
figure(12)
plot(f_g1,(abs(Y_H1)/normalize_g1))
axis([-700 700 0 0.1])
% title('Z(f) through LPF 5^t^h Order')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(12)dsb-sc_5LPF_freq','-dpng')
 
figure(13)
plot(f_g1,(abs(Y_H2)/normalize_g1))
axis([-700 700 0 0.1])
% title('Z(f) through LPF 40^t^h Order')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(13)dsb-sc_40LPF_freq','-dpng')
 
% DEMOD DSB-C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(14)
plot(t_g1,(demod_x_wc))
axis([-0.1 0.1 0 14])
% title('Demodulated Signal z_w_c(t)')
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(14)demodulated_dsb-c_time','-dpng')
 
figure(15)
plot(f_g1,abs(demod_X_wc)/normalize_g1)
axis([-700 700 0 0.1])
% title('Demodulated Signal Z_w_c(f)')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(15)demodulated_dsb-c_freq','-dpng')
 
figure(16)
plot(f_g1,(abs(Ywc_H1)/normalize_g1))
axis([-700 700 0 0.1])
% title('Z_w_c(f) through LPF 5^t^h Order')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(16)dsb-c_5LPF_freq','-dpng')
 
figure(17)
plot(f_g1,(abs(Ywc_H2)/normalize_g1))
axis([-700 700 0 0.1])
% title('Z_w_c(f) through LPF 40^t^h Order')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(17)dsb-c_40LPF_freq','-dpng')
 
% DEMODULATED AND FILTERED SIGNALS BACK IN TIME DOMAIN~~~~~~~~~~~~~~
y_H1plot = abs(ifft(Y_H1_plot));
y_H2plot = abs(ifft(Y_H2_plot));
 
ywc_H1 = ifft(Ywc_H1);
ywc_H2 = ifft(Ywc_H2);
 
figure(18)
plot(t_g1,y_H1plot-4)
axis([-0.1 0.1 -3.5 3.5])
% title('Demodulated Signal g_s_c(t) H1')
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(18)dsb-sc_LPF5_time','-dpng')
 
figure(19)
plot(t_g1,y_H2plot-4)
axis([-0.1 0.1 -3.5 3.5])
% title('Demodulated Signal g_s_c(t) H2')
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(19)dsb-sc_LPF40_time','-dpng')
 
figure(20)
plot(t_g1,abs(ywc_H1))
axis([-0.1 0.1 0 7])
%title('Demodulated Signal g_c(t) H1')
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(20)dsb-c_LPF5_time','-dpng')
 
figure(21)
plot(t_g1,abs(ywc_H2))
axis([-0.1 0.1 0 7])
% title('Demodulated Signal g_c(t) H2')
xlabel('Time (s)')
ylabel('Amplitude')
print('C:\Users\Zachary\Desktop\Matlab Plots\pt2\(21)dsb-c_LPF40_time','-dpng')
 
% END OF FILE