clc; clear variables; close all;

N = 10^6; %Number of monte carlo simulations
SNR = 0:4:40; %SNR range in dB
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
semilogy(SNR, ber1, '-r*','linewidth', 1.5); hold on;grid on;
semilogy(SNR, ber2,'-m*', 'linewidth', 1.5);hold on;grid on;
semilogy(SNR, ber3,'-c*', 'linewidth', 1.5);hold on;grid on;
%legend('User 1(FU)(BPSK) \alpha_1 = 0.8','User 2(MU)(BPSK) \alpha_2 = 0.15','User 3(NU)(BPSK) \alpha_3 = 0.05');


N = 5*10^5;

Pt = 0:4:40;			%Transmit power (dBm)
pt = (10^-3)*db2pow(Pt);	%Transmit power (linear scale)

BW = 10^6;			%Bandwidth = 1 MHz
No = -174 + 10*log10(BW);	%Noise power (dBm)			%
no = (10^-3)*db2pow(No);	%Noise power (linear scale) 

d1 = 500; d2 = 200; d3 = 70;	%Distances
a1 = 0.8; a2 = 0.15; a3 = 0.05;	%Power allocation coefficients

eta = 4;	%Path loss exponent

%Generate Rayleigh fading channel for the three users
h1 = sqrt(d1^-eta)*(randn(N/2,1) + 1i*randn(N/2,1))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(N/2,1) + 1i*randn(N/2,1))/sqrt(2);
h3 = sqrt(d3^-eta)*(randn(N/2,1) + 1i*randn(N/2,1))/sqrt(2);

%Generate noise samples for the three users
n1 = sqrt(no)*(randn(N/2,1) + 1i*randn(N/2,1))/sqrt(2);
n2 = sqrt(no)*(randn(N/2,1) + 1i*randn(N/2,1))/sqrt(2);
n3 = sqrt(no)*(randn(N/2,1) + 1i*randn(N/2,1))/sqrt(2);

%Generate random binary message data for the three users
x1 = randi([0 1],N,1);
x2 = randi([0 1],N,1);
x3 = randi([0 1],N,1);

%Create QPSKModulator and QPSKDemodulator objects
QPSKmod = comm.QPSKModulator('BitInput',true); 
QPSKdemod = comm.QPSKDemodulator('BitOutput',true); 

%Perform QPSK modulation
xmod1 = step(QPSKmod, x1);
xmod2 = step(QPSKmod, x2);
xmod3 = step(QPSKmod, x3);

%Do super position coding
x = sqrt(a1)*xmod1 + sqrt(a2)*xmod2 + sqrt(a3)*xmod3;

for u = 1:length(Pt)
	
%Received signals
   y1 = sqrt(pt(u))*x.*h1 + n1;	%At user 1
   y2 = sqrt(pt(u))*x.*h2 + n2;	%At user 2
   y3 = sqrt(pt(u))*x.*h3 + n3;	%At user 3
   
%Perform equalization
   eq1 = y1./h1;
   eq2 = y2./h2;
   eq3 = y3./h3;
   
%Decode at user 1 (Direct decoding)
   dec1 = step(QPSKdemod, eq1);
   
%Decode at user 2
   dec12 = step(QPSKdemod, eq2);		%Direct demodulation to get U1's data
   dec12_remod = step(QPSKmod, dec12);		%Remodulation of U1's data
   rem2 = eq2 - sqrt(a1*pt(u))*dec12_remod;	%SIC to remove U1's data
   dec2 = step(QPSKdemod, rem2);		%Direct demodulation of remaining signal
   
%Decode at user 3
   dec13 = step(QPSKdemod, eq3);		%Direct demodulation to get U1's data
   dec13_remod = step(QPSKmod, dec13);		%Remodulation of U1's data
   rem31 = eq3 - sqrt(a1*pt(u))*dec12_remod;	%SIC to remove U1's data
   dec23 = step(QPSKdemod, rem31);		%Direct demodulation of remaining signal to get U2's data
   dec23_remod = step(QPSKmod, dec23);		%Remodulation of U2's data
   rem3 = rem31 - sqrt(a2*pt(u))*dec23_remod;	%SIC to remove U2's data
   dec3 = step(QPSKdemod, rem3);		%Demodulate remaining signal to get U3's data
   
%BER calculation
   ber4(u) = biterr(dec1, x1)/N;
   ber5(u) = biterr(dec2, x2)/N;
   ber6(u) = biterr(dec3, x3)/N;
end

semilogy(Pt, ber4, '-k*', 'linewidth', 2); hold on; grid on;
semilogy(Pt, ber5, '-b*', 'linewidth', 2); hold on;grid on;
semilogy(Pt, ber6, '-y*', 'linewidth', 2);hold on;grid on;


%legend('User 1(QPSK) (Weakest user) \alpha_1=0.8', 'User 2(QPSK) (Middle User) \alpha_2=0.15', 'User 3 (QPSK)(Strongest user)\alpha_3=0.05');


N = 5*10^5;

Pt = 0:4:40;			%Transmit power (dBm)
pt = (10^-3)*db2pow(Pt);	%Transmit power (linear scale)

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

for u = 1:length(Pt)
	
%Received signals
   y1 = sqrt(pt(u))*x.*h1 + n1;	%At user 1
   y2 = sqrt(pt(u))*x.*h2 + n2;	%At user 2
   y3 = sqrt(pt(u))*x.*h3 + n3;	%At user 3
   
%Perform equalization
   eq1 = y1./h1;
   eq2 = y2./h2;
   eq3 = y3./h3;
   
%Decode at user 1 (Direct decoding)
   dec1 = step(QAMdemod, eq1);
   
%Decode at user 2
   dec12 = step(QAMdemod, eq2);		%Direct demodulation to get U1's data
   dec12_remod = step(QAMmod, dec12);		%Remodulation of U1's data
   rem2 = eq2 - sqrt(a1*pt(u))*dec12_remod;	%SIC to remove U1's data
   dec2 = step(QAMdemod, rem2);		%Direct demodulation of remaining signal
   
%Decode at user 3
   dec13 = step(QAMdemod, eq3);		%Direct demodulation to get U1's data
   dec13_remod = step(QAMmod, dec13);		%Remodulation of U1's data
   rem31 = eq3 - sqrt(a1*pt(u))*dec12_remod;	%SIC to remove U1's data
   dec23 = step(QAMdemod, rem31);		%Direct demodulation of remaining signal to get U2's data
   dec23_remod = step(QAMmod, dec23);		%Remodulation of U2's data
   rem3 = rem31 - sqrt(a2*pt(u))*dec23_remod;	%SIC to remove U2's data
   dec3 = step(QAMdemod, rem3);		%Demodulate remaining signal to get U3's data
   
%BER calculation
   ber7(u) = biterr(dec1, x1)/N;
   ber8(u) = biterr(dec2, x2)/N;
   ber9(u) = biterr(dec3, x3)/N;
end

semilogy(Pt, ber7, '-k+', 'linewidth', 2); hold on; grid on;
semilogy(Pt, ber8, '-g* ','linewidth', 2);
semilogy(Pt, ber9, '-r+', 'linewidth', 2);

xlabel('SNR(DB)');
ylabel('BER');
legend('User 1(FU)(BPSK) \alpha_1 = 0.8','User 2(MU)(BPSK) \alpha_2 = 0.15','User 3(NU)(BPSK) \alpha_3 = 0.05','User 1(QPSK) (Weakest user) \alpha_1=0.8', 'User 2(QPSK) (Middle User) \alpha_2=0.15', 'User 3 (QPSK)(Strongest user)\alpha_3=0.05','User 1(QAM) (Weakest user)\alpha_1 = 0.8', 'User 2(QAM)(Middle user)\alpha_2 = 0.15', 'User 3(QAM) (Strongest user)\alpha_3=0.05');

