function tx_out = tx_burst_ac(data_psdu, MSC_index, bw)
%TX_BURST Summary of this function goes here
%   Detailed explanation goes here

if bw ~= 160,
   error('Only 160MHz mode is currently supported');
end

if MSC_index<0 || MSC_index>9,
   error('Wrong MSC index');
elseif MSC_index==7 || MSC_index==9,
   error('Unsupported MSC index (supported 0...6, 8)'); 
end



%PSDU_LENGTH = 1500;
%data_psdu = randi([0 1], 1, PSDU_LENGTH*8);


% 0	BPSK	1/2
% 1	QPSK	1/2
% 2	QPSK	3/4
% 3	16-QAM	1/2
% 4	16-QAM	3/4
% 5	64-QAM	2/3
% 6	64-QAM	3/4
% 7	64-QAM	5/6 (unsupported)
% 8	256-QAM	3/4
% 9	256-QAM	5/6 (unsupported)


% system parameters (160MHz channel)
% Timing related constants
N_sd = 468;
N_sp = 16;
N_st = 484;
N_sr = 250;
N_fs = 1;
N_ss = 1;

tx_len = length(data_psdu)/8; 

if mod(length(data_psdu),8) ~= 0,
   error('PSDU length should be multiple of 8\n'); 
end

% Mode dependent constants

if MSC_index==0,
    N_bpscs = 1;
    cr = 1/2;   
    P_pattern = [1 1];  
    
elseif MSC_index==1,
    N_bpscs = 2;
    cr = 1/2;  
    P_pattern = [1 1];  
       
elseif MSC_index==2,
    N_bpscs = 2;
    cr = 3/4;
    P_pattern = [1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1];      
    
elseif MSC_index==3, 
    N_bpscs = 4;
    cr = 1/2;
    P_pattern = [1 1]; 
    
elseif MSC_index==4,
    N_bpscs = 4;
    cr = 3/4;
    P_pattern = [1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1]; 
    
elseif MSC_index==5, 
    N_bpscs = 6;
    cr = 2/3;
    P_pattern = [1 1 1 0 1 1 1 0 1 1 1 0];

elseif MSC_index==6,
    N_bpscs = 6;
    cr = 3/4;
    P_pattern = [1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1]; 

elseif MSC_index==8,
    N_bpscs = 8;
    cr = 3/4;
    P_pattern = [1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1];     
    
else
    % wrong mode
    
end

N_cbps = N_sd * N_fs * N_ss * N_bpscs;
N_dbps = N_cbps * cr;

N_DATA_SIG_B = 23;
N_TAIL_FIELD = 6;
N_PAD = 1;
N_REP = 4;
N_REP2 = 2;
N_SERVICE_FIELD = 16;
N_SCRAMBLER_STATES = 7;

trellis = poly2trellis(7,[133 171]);

% Signal symbol encoding
N_fft = 512;
N_gi = 128;
P_pilots = [-231 -203 -167 -139 -117 -89 -53 -25 25 53 89 117 139 167 203 231];
P_active = [-250:-130 -126:-6 6:126 130:250];
P_data = zeros(1, length(P_active)-length(P_pilots));
n = 1;
for k=1:length(P_active),
    if sum(P_active(k)==P_pilots) == 0,
       P_data(n) = P_active(k); % Add non-pilot carrier
       n = n + 1;
    end
end
P_data_mod = mod(P_data+N_fft, N_fft) + 1;
P_pilots_mod = mod(P_pilots+N_fft, N_fft) + 1;

LTF_left = [1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
LTF_right = [1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1];

VHTLTF80 = [LTF_left 1 LTF_right -1 -1 -1 1 1 -1 1 -1 1 1 -1 LTF_left 1 LTF_right ...
          1 -1 1 -1 0 0 0 1 -1 -1 1 ...
          LTF_left 1 LTF_right -1 -1 -1 1 1 -1 1 -1 1 1 -1 LTF_left 1 LTF_right];

VHTLTF160 = [VHTLTF80 0 0 0 0 0 0 0 0 0 0 0 VHTLTF80];      

VHTLTF = VHTLTF160;

pilot_pattern = [1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1];

% Phase rotation 
Gamma = ones(size(VHTLTF));
for k=1:length(Gamma),
   k0 = k - 251; 
   if (k0<-192),
       Gamma(k) = 1;
   elseif (k0<0 && k0>=-192),
       Gamma(k) = -1;
   elseif (k0>=0 && k0<64),
       Gamma(k) = 1; 
   elseif (k0>=64),
      Gamma(k) = -1; 
   end 
end

% Generate long training sequence
P_long = mod((-N_sr:N_sr)+N_fft, N_fft) + 1;
tx_freq = zeros(1,N_fft);
tx_freq(P_long) = VHTLTF.*Gamma;
tx_long = ifft(tx_freq);
tx_long = [tx_long(end-N_gi+1:end) tx_long tx_long(1)]; % one sample overlap

tx_long(1) = tx_long(1)*0.5;
tx_long(end) = tx_long(end)*0.5;


% Generate VHT-SIG-B
N_cbps_sig_b = N_sd * N_fs * N_ss * 1;

sig_b_bits = [0 1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 1 1 0 0 1 0]; % fixed pattenn
sig_b_bits_tail = [sig_b_bits 0 0 0 0 0 0];
sig_b_bits_all = [sig_b_bits_tail sig_b_bits_tail sig_b_bits_tail sig_b_bits_tail 0];
sig_b_bits_all = [sig_b_bits_all sig_b_bits_all];

%fprintf('%i ', sig_b_bits_all);
%fprintf('\n');

tx_sig_b_enc = convenc(sig_b_bits_all, trellis);
tx_sig_b_perm = permuter_ac(tx_sig_b_enc, N_cbps_sig_b, 1);

tx_sig_b_symb = 2*tx_sig_b_perm - 1;

tx_freq = zeros(1, N_fft);
tx_freq(P_data_mod) = tx_sig_b_symb;
tx_freq(P_pilots_mod) = pilot_pattern;

tx_freq(P_long) = tx_freq(P_long).*Gamma;

%plot(real(tx_freq(P_long)));
%hold on

tx_sig_b = ifft(tx_freq);
tx_sig_b = [tx_sig_b(end-N_gi+1:end) tx_sig_b tx_sig_b(1)]; % one sample overlap

tx_sig_b(1) = tx_sig_b(1)*0.5;
tx_sig_b(end) = tx_sig_b(end)*0.5;



% Generate data symbols
N_sym = ceil((N_SERVICE_FIELD+8*tx_len+N_TAIL_FIELD)/N_dbps);
N_data = N_sym * N_dbps;
N_pad = N_data - (N_SERVICE_FIELD+8*tx_len+N_TAIL_FIELD);

data_service = zeros(1, N_SERVICE_FIELD);
data_tail = zeros(1, N_TAIL_FIELD);
data_pad = zeros(1, N_pad);
data = [data_service data_psdu data_tail data_pad];

data_scr = scrambler(data, ones(1,N_SCRAMBLER_STATES));

data_scr(1+N_SERVICE_FIELD+tx_len*8:N_SERVICE_FIELD+tx_len*8+N_TAIL_FIELD) = data_tail;

data_enc = convenc(data_scr, trellis);

%N_cbps
%N_dbps
%N_sym

% Puncturing
data_punct = zeros(1, round(length(data_scr)/cr));
P_len = length(P_pattern);
P_nz = length(find(P_pattern == 1));
for k=1:length(data_enc)/P_len,
    data_chunk = data_enc(1+(k-1)*P_len:k*P_len);
    data_chunk = data_chunk(P_pattern==1);
    data_punct(1+(k-1)*P_nz:k*P_nz) = data_chunk;
end

% Symbol allocation, interleaving, maping, OFDM modulation
tx_data = zeros(1, N_sym*(N_fft+N_gi)+1); 
pilot_polarity = scrambler(zeros(1,N_sym+1), ones(1,N_SCRAMBLER_STATES));

for m=1:N_sym,
    
    data_sym = data_punct(1+(m-1)*N_cbps:m*N_cbps);
    data_perm = permuter_ac(data_sym, N_cbps, N_bpscs);
    tx_freq = mapper(data_perm, N_bpscs);
    tx_sym = zeros(1, N_fft);
    tx_sym(P_data_mod) = tx_freq;
    tx_sym(P_pilots_mod) = pilot_pattern * (2*pilot_polarity(m+1)-1);

    tx_sym = ifft(tx_sym);
    tx_sym = [tx_sym(end-N_gi+1:end) tx_sym tx_sym(1)]; %#ok<*AGROW> % one sample overlap

    tx_sym(1) = tx_sym(1)*0.5;
    tx_sym(end) = tx_sym(end)*0.5;
    
    tx_data(1+(m-1)*(N_fft+N_gi):m*(N_fft+N_gi)+1) = ...
        tx_data(1+(m-1)*(N_fft+N_gi):m*(N_fft+N_gi)+1) + tx_sym;   

end



N_pg = 128;

tx_out = zeros(1, (N_gi+N_fft)*(2+N_sym) + N_pg); 

% Adding Long training sequence
tx_out(1:(N_gi+N_fft)+1) = tx_out(1:(N_gi+N_fft)+1) + tx_long;

% Adding Signal symbol
tx_out(1*(N_gi+N_fft)+1:2*(N_gi+N_fft)+1) = ...
    tx_out(1*(N_gi+N_fft)+1:2*(N_gi+N_fft)+1) + tx_sig_b;

% Adding data symbols
tx_out(2*(N_gi+N_fft)+1:(2+N_sym)*(N_gi+N_fft)+1) = ...
    tx_out(2*(N_gi+N_fft)+1:(2+N_sym)*(N_gi+N_fft)+1) + tx_data;



end

