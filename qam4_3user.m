clc; clear variables; close all;

N = 5*10^5;

snr = 0:4:40;			%Transmit power (dBm)
SNR = (10^-3)*db2pow(snr);	%Transmit power (linear scale)

BW = 10^6;			%Bandwidth = 1 MHz
No = -174 + 10*log10(BW);	%Noise power (dBm)			%
no = (10^-3)*db2pow(No);	%Noise power (linear scale) 

d1 = 500; d2 = 200; d3 = 70;	%Distances
a1 = 0.8; a2 = 0.15; a3 = 0.05;	%Power allocation coefficients

eta = 4;	%Path loss exponent

%Generate Rayleigh fading channel for the three users
h1 = sqrt(d1^-eta)*(randn(N/4,1) + 1i*randn(N/4,1))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(N/4,1) + 1i*randn(N/4,1))/sqrt(2);
h3 = sqrt(d3^-eta)*(randn(N/4,1) + 1i*randn(N/4,1))/sqrt(2);

%Generate noise samples for the three users
n1 = sqrt(no)*(randn(N/4,1) + 1i*randn(N/4,1))/sqrt(2);
n2 = sqrt(no)*(randn(N/4,1) + 1i*randn(N/4,1))/sqrt(2);
n3 = sqrt(no)*(randn(N/4,1) + 1i*randn(N/4,1))/sqrt(2);

%Generate random binary message data for the three users
x1 = randi([0 1],N/2,1);
x2 = randi([0 1],N/2,1);
x3 = randi([0 1],N/2,1);

%Create QAMModulator and QAMModulator objects
QAMmod = comm.RectangularQAMModulator(4,'BitInput',true); 
QAMdemod = comm.RectangularQAMDemodulator(4,'BitOutput',true);
%Perform QPSK modulation
xmod1 = step(QAMmod, x1);
xmod2 = step(QAMmod, x2);
xmod3 = step(QAMmod, x3);

%Do super position coding
x = sqrt(a1)*xmod1 + sqrt(a2)*xmod2 + sqrt(a3)*xmod3;

for u = 1:length(snr)
	
%Received signals
   y1 = sqrt(SNR(u))*x.*h1 + n1;	%At user 1
   y2 = sqrt(SNR(u))*x.*h2 + n2;	%At user 2
   y3 = sqrt(SNR(u))*x.*h3 + n3;	%At user 3
   
%Perform equalization
   eq1 = y1./h1;
   eq2 = y2./h2;
   eq3 = y3./h3;
   
%Decode at user 1 (Direct decoding)
   dec1 = step(QAMdemod, eq1);
   
%Decode at user 2
   dec12 = step(QAMdemod, eq2);		%Direct demodulation to get U1's data
   dec12_remod = step(QAMmod, dec12);		%Remodulation of U1's data
   rem2 = eq2 - sqrt(a1*SNR(u))*dec12_remod;	%SIC to remove U1's data
   dec2 = step(QAMdemod, rem2);		%Direct demodulation of remaining signal
   
%Decode at user 3
   dec13 = step(QAMdemod, eq3);		%Direct demodulation to get U1's data
   dec13_remod = step(QAMmod, dec13);		%Remodulation of U1's data
   rem31 = eq3 - sqrt(a1*SNR(u))*dec12_remod;	%SIC to remove U1's data
   dec23 = step(QAMdemod, rem31);		%Direct demodulation of remaining signal to get U2's data
   dec23_remod = step(QAMmod, dec23);		%Remodulation of U2's data
   rem3 = rem31 - sqrt(a2*SNR(u))*dec23_remod;	%SIC to remove U2's data
   dec3 = step(QAMdemod, rem3);		%Demodulate remaining signal to get U3's data
   
%BER calculation
   ber1(u) = biterr(dec1, x1)/N;
   ber2(u) = biterr(dec2, x2)/N;
   ber3(u) = biterr(dec3, x3)/N;
end

semilogy(snr, ber1, '-r*', 'linewidth', 2); hold on; grid on;
semilogy(snr, ber2, '-m* ','linewidth', 2);
semilogy(snr, ber3, '-b*', 'linewidth', 2);

xlabel('SNR (dB)');
ylabel('BER');
legend('User 1 (Weakest user)\alpha_1 = 0.8', 'User 2(Middle user)\alpha_2 = 0.15', 'User 3 (Strongest user)\alpha_3=0.05');
title('4-QAM for 3-USER');