%BPSK simulation using a carrier cosine wave with ISI
clc;
close all;
clear all;
%figure(1)
n=160
for i=1:n
    data(i)= 2*round(rand)-1;
end
%create modulated BPSK signal
%first expand the bit stream
exdata=[];
for i=1:length(data)
  for rep=1:5
  exdata= [exdata data(i)];
  end
end
ts=.1;
t=1:ts:80.9;
carrier=cos(pi*t);
%multiply expanded bitstream by cosine wave with carrier frequency
%this is the BPSK that is to be transmitted over the channel
bpsk=carrier.*exdata;
%bpsk=[bpsk(length(bpsk)-1) bpsk(length(bpsk)) bpsk];
%plot(bpsk)
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
for k=5:5:795
    for l = 1:5
    %al=round(rand)+j*round(rand)
     rec(k+l)=bpsk(k+l)+al*bpsk(k-5+l);
    end
end
rxdata=rec+ no ;
%begin demodulation
%first multiply recieved bitstream by cosine wave with carrier frequency
%figure(2)
uncarry=rxdata.*carrier;
%plot(uncarry)
%demodulate by integrating 
dec1=[];
for inc=1:5:length(uncarry)   
  dec=trapz(inc:inc+4,uncarry(inc:inc+4));
  dec1=[dec1 dec];
end
%make decision with a threshold of zero
demod=[];
for i=1:length(dec1)
    if dec1(i)>0
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
error
ber=error/n
figure(3)
title('Comparing the Bits at transmitter and receiver')
stem(data)
hold
stem(demod,'rx')