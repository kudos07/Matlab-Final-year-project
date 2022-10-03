
N = 5*10^5;

pt = 0:2:40;	
Pt = (10^-3)*db2pow(pt);	

BW = 10^6;	
No = -174 + 10*log10(BW);	
no = (10^-3)*db2pow(No);

x1 = randi([0 1],1,N); 
x2 = randi([0 1],1,N); 

xmod1 = 2*x1 - 1;
xmod2 = 2*x2 - 1;

%Set power weights for users
a1 = 0.8; a2 = 0.2;

%Do superposition coding
x = sqrt(a1)*xmod1 + sqrt(a2)*xmod2 ;

for u = 1:length(Pt)
y1 = awgn(x,pt(u),'measured'); 
y2 = awgn(x,pt(u),'measured'); 

x1_hat = ones(1,N);
x1_hat(y1 < 0) = 0;

% at user 2
x11_est = ones(1,N); 
x11_est(y2 < 0) = 0;
x11_est(x11_est == 0) = -1;                                                               

rem2 = y2 - sqrt(a1)*x11_est;

x2_hat = zeros(1,N);
x2_hat(rem2>0) = 1;

%Estimate BER
ber1(u) = biterr(x1,x1_hat)/N;
ber2(u) = biterr(x2,x2_hat)/N;
end

%plot BER curves
semilogy(pt, ber1,'-k*', 'linewidth', 2);hold on;grid on;
semilogy(pt, ber2, '-c*', 'linewidth', 2);hold on;grid on;

N = 5*10^5;

pt = 0:4:40;			%Tsnr (dBm)
Pt = (10^-3)*db2pow(pt);	% (linear scale)

BW = 10^6;			%Bandwidth = 1 MHz
No = -174 + 10*log10(BW);	%Noise power (dBm)
no = (10^-3)*db2pow(No);	%Noise power (linear scale) 

d1 = 500; d2 = 200;	%Distances
a1 = 0.9; a2=0.1;	%Power allocation coefficients

eta = 4;	%Path loss exponent

%Generate Rayleigh fading channel for the three users
h1 = sqrt(d1^-eta)*(randn(N/2,1) + 1i*randn(N/2,1))/sqrt(2) ;
h2 = sqrt(d2^-eta)*(randn(N/2,1) + 1i*randn(N/2,1))/sqrt(2);


%Generate noise samples for the three users
n1 = sqrt(no)*(randn(N/2,1) + 1i*randn(N/2,1))/sqrt(2);
n2 = sqrt(no)*(randn(N/2,1) + 1i*randn(N/2,1))/sqrt(2);


%Generate random binary message data for the three users
x1 = randi([0 1],N,1);
x2 = randi([0 1],N,1);


%Create QPSKModulator and QPSKDemodulator objects
QPSKmod = comm.QPSKModulator('BitInput',true); 
QPSKdemod = comm.QPSKDemodulator('BitOutput',true); 

%Perform QPSK modulation
xmod1 = step(QPSKmod, x1);
xmod2 = step(QPSKmod, x2);


%Do super position coding
x = sqrt(a1)*xmod1 + sqrt(a2)*xmod2;

for u = 1:length(pt)
	
%Received signals
   y1 = sqrt(Pt(u))*x.*h1 + n1;	%At user 1
   y2 = sqrt(Pt(u))*x.*h2 + n2;	%At user 2
   
%Perform equalization
   eq1 = y1./h1;
   eq2 = y2./h2;
   
%Decode at user 1 (Direct decoding)
   dec1 = step(QPSKdemod, eq1);
   
%Decode at user 2
   dec12 = step(QPSKdemod, eq2);		%Direct demodulation to get U1's data
   dec12_remod = step(QPSKmod, dec12);		%Remodulation of U1's data
   rem2 = eq2 - sqrt(a1*Pt(u))*dec12_remod;	%SIC to remove U1's data
   dec2 = step(QPSKdemod, rem2);		%Direct demodulation of remaining signal
   
   
   
%BER calculation
   ber3(u) = biterr(dec1, x1)/N;
   ber4(u) = biterr(dec2, x2)/N;
  
end

semilogy(pt, ber3,'-g*', 'linewidth', 2); hold on;grid on;

semilogy(pt, ber4, '-b*', 'linewidth', 2);hold on;grid on;





N = 5*10^5;

Pt = 0:4:40;			%Transmit power (dBm)
pt = (10^-3)*db2pow(Pt);	%Transmit power (linear scale)

BW = 10^6;			%Bandwidth = 1 MHz
No = -174 + 10*log10(BW);	%Noise power (dBm)			%
no = (10^-3)*db2pow(No);	%Noise power (linear scale) 

d1 = 500; d2 = 200; 
a1 = 0.91; a2 = 0.09; 

eta = 4;	%Path loss exponent

%Generate Rayleigh fading channel for the three users
h1 = sqrt(d1^-eta)*(randn(N/4,1) + 1i*randn(N/4,1))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(N/4,1) + 1i*randn(N/4,1))/sqrt(2);

%Generate noise samples for the three users
n1 = sqrt(no)*(randn(N/4,1) + 1i*randn(N/4,1))/sqrt(2);
n2 = sqrt(no)*(randn(N/4,1) + 1i*randn(N/4,1))/sqrt(2);


%Generate random binary message data for the three users
x1 = randi([0 1],N/2,1);
x2 = randi([0 1],N/2,1);


%Create QAMModulator and QAMDemodulator objects
QAMmod = comm.RectangularQAMModulator(4,'BitInput',true); 
QAMdemod = comm.RectangularQAMDemodulator(4,'BitOutput',true);

%Perform QAM modulation
xmod1 = step(QAMmod, x1);
xmod2 = step(QAMmod, x2);


%Do super position coding
x = sqrt(a1)*xmod1 + sqrt(a2)*xmod2;

for u = 1:length(Pt)
	
%Received signals
   y1 = sqrt(pt(u))*x.*h1 + n1;	%At user 1
   y2 = sqrt(pt(u))*x.*h2 + n2;	%At user 2
 
   
%Perform equalization
   eq1 = y1./h1;
   eq2 = y2./h2;
  
   
%Decode at user 1 (Direct decoding)
   dec1 = step(QAMdemod, eq1);
   
%Decode at user 2
   dec12 = step(QAMdemod, eq2);		%Direct demodulation to get U1's data
   dec12_remod = step(QAMmod, dec12);		%Remodulation of U1's data
   rem2 = eq2 - sqrt(a1*pt(u))*dec12_remod;	%SIC to remove U1's data
   dec2 = step(QAMdemod, rem2);		%Direct demodulation of remaining signal
   
%BER calculation
   ber5(u) = biterr(dec1, x1)/N;
   ber6(u) = biterr(dec2, x2)/N;
  
end

semilogy(Pt, ber5,'-m*' ,'linewidth', 2); hold on; grid on;
semilogy(Pt, ber6, '-r*', 'linewidth', 2);hold on;grid on;
legend( 'User 1(BPSK) (FU) \alpha_1 = 0.8','User 2(BPSK) (NU) \alpha_2 = 0.2','User 1(QPSK) (FU) \alpha_1 = 0.91','User 2(QPSK) (NU) \alpha_2 = 0.09','User 1(FU)(QAM) \alpha_1 = 0.9','User 2(NU)(QAM) \alpha_2 = 0.1');
xlabel('SNR (dB)');
ylabel('BER');

