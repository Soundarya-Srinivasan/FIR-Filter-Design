clc;
clear all;
close all;

fs = 8000;
f = [0 600 1000 1400 4000 8000];
a = [1.25 2 1 .25 0.1 0];
order = 31;
bits = 16;
intgr = 9;
A = 2;

coeffs = fir2(order, f/fs, a);
quantizedCoeffs = sfi(coeffs,bits,bits-intgr-1);

ts = 1/fs;
ns = 256;
t = 0:ts:ts*(ns-1);

x = chirp(t,0,ts*ns,fs/2);

% Non Qauntized
figure;
X = abs(fft(x,ns));
X = X(1:length(X)/2);
frq = 1:1:length(X);

subplot(3,1,1);
plot(frq*(fs/ns),X);
title("Non Quantized");
grid on;
%plot normalized frequency of filter

[h,w] = freqz(coeffs,A,256);
subplot(3,1,2);
plot(w/(2*pi),10*log(abs(h)));
grid on;
%plot fft of filtered signal
Y = filter(coeffs,A,x);
Y = abs(fft(Y,ns));
Y = Y(1:length(Y)/2);

frq = 1:1:length(Y);
subplot(3,1,3);
plot(frq*(fs/ns),Y);
grid on;

%Quantized 
figure;
qX = abs(fft(x,ns));
qX = qX(1:length(qX)/2);
qfrq = 1:1:length(qX);

subplot(3,1,1);
plot(qfrq*(fs/ns),qX);
title("Quantized");
grid on;

[qh,qw] = freqz(quantizedCoeffs.data,A,256);
subplot(3,1,2);
plot(qw/(2*pi),10*log(abs(qh)));

grid on;

%plot fft of filtered signal
qY = filter(quantizedCoeffs.data,A,x);
qY = abs(fft(qY,ns));
qY = qY(1:length(qY)/2);

qfrq = 1:1:length(qY);
subplot(3,1,3);
plot(qfrq*(fs/ns),qY);
grid on;

error = immse(Y, qY);
quantizedData = quantizedCoeffs.data; 
minMaxNormCoeff = (quantizedData-min(quantizedData))./(max(quantizedData)-min(quantizedData));
NormCoeff = ((2^15 )- 1)*minMaxNormCoeff;
