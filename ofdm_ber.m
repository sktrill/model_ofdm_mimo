% Vinay Mohan Jonnakuti
% Wireless Communication:
% BPSK modulation and demodulation with increased bit rate
%BPSK simulation using a carrier cosine wave with increased bit rate
clc;
close all;
clear all;

% no of sub carrier channels
channels=6;
c=6;
% bits per channels
bits=54
% total no of bits to be transmitted
n=324;
for i=1:n
    data(i)= 2*round(rand)-1;
end

% Converting the series into parallel for the channels
s = reshape(data,c,bits);
% BPSK modulation for a single channel
signal1 = s(1,:);
%first expand the bit stream
exdata1=[];
exdata2=[];
exdata3=[];
exdata4=[];
exdata5=[];
exdata6=[];
for i=1:length(signal1)
  for rep=1:2
  exdata1= [exdata1 signal1(i)];
  end
end

signal2 = s(2,:);
%first expand the bit stream
exdata=[];
for i=1:length(signal2)
  for rep=1:2
  exdata2= [exdata2 signal2(i)];
  end
end

signal3= s(3,:);   
%first expand the bit stream
exdata=[];
for i=1:length(signal3)
  for rep=1:2
  exdata3= [exdata3 signal3(i)];
  end
end

signal4= s(4,:);
%first expand the bit stream
exdata=[];
for i=1:length(signal4)
  for rep=1:2
  exdata4= [exdata4 signal4(i)];
  end
end

signal5= s(5,:);
%first expand the bit stream
exdata=[];
for i=1:length(signal5)
  for rep=1:2
  exdata5= [exdata5 signal5(i)];
  end
end

signal6= s(6,:);
%first expand the bit stream
exdata=[];
for i=1:length(signal6)
  for rep=1:2
  exdata6= [exdata6 signal6(i)];
  end
end

% Bpsk modulation
m=10*n;
% Generating the carrier signal  
ts=.1;
tp=1:ts:11.79;
carrier1=cos(2*pi*tp);
% Generating the modulated signal 1
bpsk_sig1=exdata1.*carrier1; 
% Generating the modulated signal 2 
carrier2=cos(4*pi*tp);
bpsk_sig2=exdata2.*carrier2; 
% Generating the modulated signal 3 
carrier3=cos(6*pi*tp);
bpsk_sig3=exdata3.*carrier3;
% Generating the modulated signal 4
carrier4=cos(8*pi*tp);
bpsk_sig4=exdata4.*carrier4;
% Generating the modulated signal 5
carrier5=cos(10*pi*tp);
bpsk_sig5=exdata5.*carrier5;
carrier6=cos(12*pi*tp);
% Generating the modulated signal
bpsk_sig6=exdata6.*carrier6;
% taking the iFFT of each of these signals
if_sig1=ifft(bpsk_sig1);
if_sig2=ifft(bpsk_sig2);
if_sig3=ifft(bpsk_sig3);
if_sig4=ifft(bpsk_sig4);
if_sig5=ifft(bpsk_sig5);
if_sig6=ifft(bpsk_sig6);

fin(1,:)=if_sig1;
fin(2,:)=if_sig2;
fin(3,:)=if_sig3;
fin(4,:)=if_sig4;
fin(5,:)=if_sig5;
fin(6,:)=if_sig6;

transmit=reshape(fin,1,648);


% generating the noise
% p=rand(1,800)*2*pi;
p=rand*2*pi;
snr=10;
r=sqrt(-1*(1/snr*log(1 - rand)));
% no = 5*(r.* exp(j*p));
no = (r.* exp(j*p));    
% value of alpha 
al=rand+j*rand;
%al=1;
% Spreading channel with the alpha as the variable

for k=2:2:646
    for l = 1:2
    %al=round(rand)+j*round(rand)
     rec(k+l)=transmit(k+l)+al*transmit(k-2+l);
    end
end

rxdata=rec+ no ;

% Converting from serial to parallel
myrec=reshape(rxdata,6,108)

rxdata1=fft(myrec(1,:));
rxdata2=fft(myrec(2,:));
rxdata3=fft(myrec(3,:));
rxdata4=fft(myrec(4,:));
rxdata5=fft(myrec(5,:));
rxdata6=fft(myrec(6,:));

% taking the FFT


%begin demodulation
%first multiply recieved bitstream by cosine wave with carrier frequency


uncarry1=rxdata1.*carrier1;
uncarry2=rxdata2.*carrier2;
uncarry3=rxdata3.*carrier3;
uncarry4=rxdata4.*carrier4;
uncarry5=rxdata5.*carrier5;
uncarry6=rxdata6.*carrier6;

%plot(uncarry)
%demodulate by integrating 
dec1=[];
dec2=[];
dec3=[];
dec4=[];
dec5=[];
dec6=[];
for inc=1:2:length(uncarry1)   
  dec=trapz(inc:inc+1,uncarry1(inc:inc+1));
  dec1=[dec1 dec];
end
%2
for inc=1:2:length(uncarry2)   
  dec=trapz(inc:inc+1,uncarry2(inc:inc+1));
  dec2=[dec2 dec];
end
%3
for inc=1:2:length(uncarry3)   
  dec=trapz(inc:inc+1,uncarry3(inc:inc+1));
  dec3=[dec3 dec];
end
%4
for inc=1:2:length(uncarry4)   
  dec=trapz(inc:inc+1,uncarry4(inc:inc+1));
  dec4=[dec4 dec];
end
%5
for inc=1:2:length(uncarry5)   
  dec=trapz(inc:inc+1,uncarry5(inc:inc+1));
  dec5=[dec5 dec];
end
%6
for inc=1:2:length(uncarry6)   
  dec=trapz(inc:inc+1,uncarry6(inc:inc+1));
  dec6=[dec6 dec];
end
final_rec(1,:)=dec1;
final_rec(2,:)=dec2;
final_rec(3,:)=dec3;
final_rec(4,:)=dec4;
final_rec(5,:)=dec5;
final_rec(6,:)=dec6;

fin_rec_parallel=reshape(final_rec,1,324)
%make decision with a threshold of zero
demod=[];
for i=1:length(fin_rec_parallel)
    if fin_rec_parallel(i)>0
        demod=[demod 1];
    else
        demod=[demod -1];
    end
end

%stem(demod)

%calculate errors
error=0;
for i=1:length(demod)
    if data(i)~=demod(i)
        error=error+1;
    end
end
error;

ber=error/324
figure(3)
stem(data)
hold
stem(demod,'rx')