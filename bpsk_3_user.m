clc; clear variables; close all;

N = 10^6; %Number of monte carlo simulations
SNR = 0:2:40; %SNR range in dB
snr = db2pow(SNR); %SNR range in linear scale

%Generate random data bits for transmission
x1 = randi([0 1],1,N); %Data bits of user 1
x2 = randi([0 1],1,N); %Data bits of user 2
x3 = randi([0 1],1,N); %Data bits of user 3

%Do BPSK modulation of data
xmod1 = 2*x1 - 1;
xmod2 = 2*x2 - 1;
xmod3 = 2*x3 - 1;

%Set power weights for users
a1 = 0.8; a2 = 0.15; a3 = 0.05;

%Do superposition coding
x = sqrt(a1)*xmod1 + sqrt(a2)*xmod2 + sqrt(a3)*xmod3;

%Add AWGN to x (Transmit x through an AWGN channel)
for u = 1:length(snr)
y1 = awgn(x,SNR(u),'measured'); %Received signal at user 1 corrupted by AWGN
y2 = awgn(x,SNR(u),'measured'); %Received signal at user 2 corrupted by AWGN
y3 = awgn(x,SNR(u),'measured'); %Received signal at user 3 corrupted by AWGN

%AT USER 1
%Direct decoding of x from y1
x1_hat = ones(1,N); %just a buffer
x1_hat(y1 < 0) = 0;

%AT USER 2
%Direct decoding of x from y2
x11_est = ones(1,N); %just a buffer
x11_est(y2 < 0) = 0; %Estimate user 1's signal first
x11_est(x11_est == 0) = -1; %Remodulate user 1's signal
%Subtract remodulated x11_est component from y2
rem2 = y2 - sqrt(a1)*x11_est;

%Decode x2 from rem
x2_hat = zeros(1,N);
x2_hat(rem2>0) = 1;



%AT USER 3
%Direct decoding of x from y3
x11_est2 = ones(1,N); %just a buffer
x11_est2(y3 < 0) = 0; %Estimate user 1's signal first
x11_est2(x11_est2 == 0) = -1; %Remodulate user 1's signal
%Subtract remodulated x11_est2 component from y3
rem22 = y3 - sqrt(a1)*x11_est2;


x11_est3 = ones(1,N); %just a buffer
x11_est3(rem22 < 0) = 0; %Estimate user 2's signal first
x11_est3(x11_est3 == 0) = -1; %Remodulate user 2's signal
%Subtract remodulated x11_est3 component from rem22
rem3 = rem22 - sqrt(a2)*x11_est3;

%Decode x2 from rem
x3_hat = zeros(1,N);
x3_hat(rem3>0) = 1;

%Estimate BER
ber1(u) = biterr(x1,x1_hat)/N;
ber2(u) = biterr(x2,x2_hat)/N;
ber3(u) = biterr(x3,x3_hat)/N;
end

%plot BER curves
semilogy(SNR, ber1, '-r*','linewidth', 1.5); hold on;
semilogy(SNR, ber2,'-m*', 'linewidth', 1.5);
semilogy(SNR, ber3,'-c*', 'linewidth', 1.5);grid on;
legend('User 1(FU) \alpha_1 = 0.8','User 2(MU) \alpha_2 = 0.15','User 3(NU) \alpha_3 = 0.05');
xlabel('SNR (dB)');
ylabel('BER');
title('BPSK 3 USER');