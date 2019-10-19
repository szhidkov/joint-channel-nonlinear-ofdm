function tx_out = multipath_channel_ac(tx_in, model, snr, L)
%MULTIPATH_CHANNEL Summary of this function goes here
%   Detailed explanation goes here

dt = 6.25 / L;

% AWGN:
% Single path channel

% Model A: 
% Corresponds to a typical office environment 
% for NLOS conditions and 50ns average rms delay spread

% Model B:
% Corresponds to a typical large open space and 
% office environments for NLOS conditions and 
% 100ns average rms delay spread

% Model C: 
% Corresponds to a typical large open space environment 
% for NLOS conditions and 150ns average rms delay spread

if strcmp(model, 'AWGN'),
    del = 0;
    att = 0;

elseif strcmp(model, 'A'),

    del = [0 10 20 30 40 50 60 70 80 90 110 140 170 200 240 290 340 390]; % ns
    att = [0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 -4.7 -7.3 ...
           -9.9 -12.5 -13.7 -18.0 -22.4 -26.7]; % dB

elseif strcmp(model, 'B'),
    
    del = [0 10 20 30 50 80 110 140 180 230 280 330 380 430 490 560 640 730];
    att = [-2.6 -3.0 -3.5 -3.9 0.0 -1.3 -2.6 -3.9 -3.4 -5.6 -7.7 -9.9 ...
           -12.1 -14.3 -15.4 -18.4 -20.7 -24.6];

elseif strcmp(model, 'C'),       
    del = [0 10 20 30 50 80 110 140 180 230 280 330 400 490 600 730 880 1050]; 
    att = [-3.3 -3.6 -3.9 -4.2 0.0 -0.9 -1.7 -2.6 -1.5 -3.0 -4.4 ...
        -5.9 -5.3 -7.9 -9.4 -13.2 -16.3 -21.2];

end

% Generate random channel realization
N_taps = round(max(del)/dt) + 1;
h = zeros(1,N_taps);
for k=1:length(del),
   t = 1 + round(del(k)/dt);
   a = 10^(att(k)*0.05);
   hk = a*(randn + 1i*randn);
   h(t) = h(t) + hk; 
end

% plot(abs(h))

tx_out = filter(h,1,tx_in);
tx_out = tx_out / std(tx_out);
%tx_out = tx_in;

snr_eq = snr - 10*log10(L);

Sw = 10^(-snr_eq*0.05);
tx_out = tx_out + 1/sqrt(2)*Sw*(randn(size(tx_out)) + 1i*randn(size(tx_out)));

end

