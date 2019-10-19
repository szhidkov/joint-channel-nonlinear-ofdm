function [d3, d5] = calcd(X, N_fft)
%CALCD Summary of this function goes here
%   Detailed explanation goes here

S = X(1+mod(N_fft/2+(1:N_fft),N_fft));
S_conj = conj(S(end:-1:1)); 
d3 = conv(S, S);
d3 = conv(d3, S_conj);
d5 = conv(d3, S);
d5 = conv(d5, S_conj);
d3 = d3((N_fft-1)+1:2*N_fft-1);
d5 = d5((2*N_fft-2)+1:3*N_fft-2);
d3a = zeros(1,N_fft);
d3a(1+mod(N_fft/2+(1:N_fft),N_fft)) = d3;
d5a = zeros(1,N_fft);
d5a(1+mod(N_fft/2+(1:N_fft),N_fft)) = d5;
d3 = d3a;
d5 = d5a;


end

