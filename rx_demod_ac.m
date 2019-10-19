function [data_out, c_hat] = rx_demod_ac(rx_in, MSC_index, bw, PSDU_LENGTH, enable_ndc, c_hat, n_iter)
%RX_BURST1 Summary of this function goes here
%   Detailed explanation goes here

if bw ~= 160,
   error('Only 160MHz mode is currently supported');
end

if MSC_index<0 || MSC_index>9,
   error('Wrong MSC index');
elseif MSC_index==7 || MSC_index==9,
   error('Unsupported MSC index (supported 0...6, 8)'); 
end

if n_iter==1,
    alp = 1.0;
elseif n_iter==2,
    alp = [1.0 1.0];
elseif n_iter==3,
    alp = [0.8 1.0 1.0];
else
    alp = ones(1,n_iter);
    alp(1) = 0.5;
end
                
% Mode dependent settings
if MSC_index==6,
    lambda = 0;
    gamma_sync = 0.5;
    gamma_data = 0.25;  % 0.15 best with fec G_fec=0  / %0.25 best with pure slicer
    N_itr = 15;
    R_th = 0.004;
elseif MSC_index==8,
    lambda = 0;
    gamma_sync = 0.5;
    gamma_data = 0.25;  % 0.15 best with fec G_fec=0  / %0.25 best with pure slicer
    N_itr = 15;
    R_th = 0.004;
else
    lambda = 0;
    gamma_sync = 0.5;
    gamma_data = 0.2;
    N_itr = 15;  
    R_th = 0.004;
end

trellis = poly2trellis(7,[133 171]);

b_rx = [-4.749621e-003 3.645163e-003 5.962925e-003 -3.269992e-003 -7.277667e-003 2.640417e-003 8.674519e-003 -1.713875e-003 -1.013080e-002 4.424068e-004 1.162100e-002 1.230960e-003 -1.311742e-002 -3.378575e-003 1.459091e-002 6.099455e-003 -1.601163e-002 -9.539999e-003 1.734992e-002 1.393414e-002 -1.857709e-002 -1.968748e-002 1.966632e-002 2.757358e-002 -2.059338e-002 -3.926643e-002 2.133736e-002 5.915005e-002 -2.188133e-002 -1.033642e-001 2.221281e-002 3.173915e-001 4.776758e-001 3.173915e-001 2.221281e-002 -1.033642e-001 -2.188133e-002 5.915005e-002 2.133736e-002 -3.926643e-002 -2.059338e-002 2.757358e-002 1.966632e-002 -1.968748e-002 -1.857709e-002 1.393414e-002 1.734992e-002 -9.539999e-003 -1.601163e-002 6.099455e-003 1.459091e-002 -3.378575e-003 -1.311742e-002 1.230960e-003 1.162100e-002 4.424068e-004 -1.013080e-002 -1.713875e-003 8.674519e-003 2.640417e-003 -7.277667e-003 -3.269992e-003 5.962925e-003 3.645163e-003 -4.749621e-003];
b_rx_cplx = b_rx+1i*b_rx;
b_delay = (length(b_rx) - 1)/2;

rx_flt = filter(b_rx_cplx, 1, rx_in);

rx_dec = rx_flt(b_delay+1:2:end);

%rx_dec = tx_out;

% Channel estimation (simple LS + smoothing or IFFT interp)
% 0	BPSK	1/2
% 1	QPSK	1/2
% 2	QPSK	3/4
% 3	16-QAM	1/2
% 4	16-QAM	3/4
% 5	64-QAM	2/3
% 6	64-QAM	3/4
% 7	64-QAM	5/6
% 8	256-QAM	3/4
% 9	256-QAM	5/6
% MSC_index = 6;

% system parameters (160MHz channel)
% Timing related constants
N_sd = 468;
N_sp = 16;
N_st = 484;
N_sr = 250;
N_fs = 1;
N_ss = 1;

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

if N_bpscs==2,
    scl = 1/sqrt(2);
elseif N_bpscs==4,
    scl = 1/sqrt(10);
elseif N_bpscs==6,
    scl = 1/sqrt(42);
elseif N_bpscs==8,
    scl = 1/sqrt(170);    
else
    scl = 1;
end

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
P_active_mod = mod(P_active+N_fft, N_fft) + 1;

P_long = mod((-N_sr:N_sr)+N_fft, N_fft) + 1;

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

N_long = length(P_active);

T3 = (2*N_long-1)*2/3*(2^N_bpscs-1)*scl^2;
T5 = 6*N_long*(N_long-1)*(2/3*(2^N_bpscs-1)*scl^2)^2 - (3*N_long-4)*(4/45*(2^N_bpscs-1)*(7*2^N_bpscs-13))*scl^4; 
T3_bpsk = (2*N_long-1);
T5_bpsk = 6*N_long*(N_long-1) - (3*N_long-4); 

N_sym = ceil((N_SERVICE_FIELD+8*PSDU_LENGTH+N_TAIL_FIELD)/N_dbps);
N_data = N_sym * N_dbps;
N_pad = N_data - (N_SERVICE_FIELD+8*PSDU_LENGTH+N_TAIL_FIELD);

F1 = fft(rx_dec(N_gi+1:N_gi+N_fft));


R_ref = zeros(1,N_fft);
R_ref(P_long) = VHTLTF.*Gamma;



H_est1 = zeros(1,N_fft);
H_est1(P_long) = F1(P_long)./R_ref(P_long);

H = H_est1;

%plot(H,'r.');
%stop

% Equalization and detection for SIGNAL symbol
R_sig = fft(rx_dec(2*N_gi+N_fft+1:2*N_gi+2*N_fft));
R_sig_eq = zeros(size(R_sig));
R_sig_eq(P_active_mod) = R_sig(P_active_mod) ./ H(P_active_mod);

R_sig_eq(P_long) = R_sig_eq(P_long) .* Gamma;

%plot(real(R_sig_eq(P_long)) ,'r:');
%stop

data_rx = R_sig_eq(P_data_mod);
csi_rx = H(P_data_mod);


% SIGNAL field demodulation/decoding
N_cbps_sig_b = N_sd * N_fs * N_ss * 1;
llr = demapper_csi(data_rx, csi_rx, 1);
llr_out = depermuter_ac(llr, N_cbps_sig_b, 1);
dec = vitdec(llr_out, trellis, 21, 'term', 'unquant');

% 
dec = [0 1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0];



% Reconstruction of SIGNAL symbol (in frequency domain)
tx_sig_b_enc = convenc(dec, trellis);
tx_sig_b_perm = permuter_ac(tx_sig_b_enc, N_cbps_sig_b, 1);
tx_sig_b_symb = 2*tx_sig_b_perm - 1;
tx_freq = zeros(1, N_fft);
tx_freq(P_data_mod) = tx_sig_b_symb;
tx_freq(P_pilots_mod) = pilot_pattern;
tx_freq(P_long) = tx_freq(P_long).*Gamma;



% Data-aided Channel estimation
H_est2 = zeros(1,N_fft);
H_est2(P_long) = R_sig(P_long)./tx_freq(P_long);
H = 1/2*H + 1/2*H_est2;


pilot_polarity = scrambler(zeros(1,N_sym+1), ones(1,N_SCRAMBLER_STATES));



if (enable_ndc),
    
    % nonlinear distortion compensation
    X2_hat = tx_freq;
    [d3, d5] = calcd(X2_hat, N_fft);
    U2 = [d3(P_active_mod) - T3*X2_hat(P_active_mod); d5(P_active_mod) - T5*X2_hat(P_active_mod)].';
    
    X1_hat = R_ref;
    [d3, d5] = calcd(X1_hat, N_fft);
    U1 = [d3(P_active_mod) - T3_bpsk*X1_hat(P_active_mod); d5(P_active_mod) - T5_bpsk*X1_hat(P_active_mod)].';   
    
    F_lts = F1;
    
    var4 = zeros(1,2*N_itr+1);
    
    var4(1) = std(F1(P_active_mod) - H(P_active_mod).*X1_hat(P_active_mod))^2 ...
            + std(R_sig(P_active_mod) - H(P_active_mod).*X2_hat(P_active_mod))^2;
    
    c_hat_old = zeros(2,1); 
    H_old = H;
        
    for itr=1:N_itr,
        
           
        Am = diag([H(P_active_mod) H(P_active_mod)])' * ...
             diag([H(P_active_mod) H(P_active_mod)]);
        Um = [U1; U2];
        Rm = [F1(P_active_mod) R_sig(P_active_mod)];
        Hm = [H(P_active_mod) H(P_active_mod)];
        Xm = [X1_hat(P_active_mod) X2_hat(P_active_mod)];
        
        c_hat3 = inv(Um'*Am*Um) * Um' * Am * (Rm./Hm - Xm).';
        c_hat_new = c_hat3;   
        
        
        var4_new = std( F1(P_active_mod) - H(P_active_mod).*(X1_hat(P_active_mod) +(U1*c_hat_new).') )^2 + ...
              std( R_sig(P_active_mod) - H(P_active_mod).*(X2_hat(P_active_mod) +(U2*c_hat_new).') )^2;
        
        if var4_new > var4(2*itr-1),
           c_hat_new = c_hat_old;
           if itr<3,
              fprintf('%i', itr*2-1);
           end
           break;
        else
            c_hat_old = c_hat_new;
        end
        
        var4(2*itr) = var4_new;
        
        H_new = zeros(1, N_fft);
        H_new(P_active_mod) =    1/2 * F1(P_active_mod)./(X1_hat(P_active_mod)+(U1*c_hat_new).') + ...
                           1/2 * R_sig(P_active_mod)./(X2_hat(P_active_mod)+(U2*c_hat_new).'); 
        
        var4_new = ...
              std( F1(P_active_mod) - H_new(P_active_mod).*(X1_hat(P_active_mod) +(U1*c_hat_new).') )^2 + ...
              std( R_sig(P_active_mod) - H_new(P_active_mod).*(X2_hat(P_active_mod) +(U2*c_hat_new).') )^2;
        
        if var4_new > var4(2*itr),
            H = H_old;
            if itr<3,
                fprintf('%i', itr*2);
            end
            break;
        else
            H_old = H_new;
        end
        
        var4(2*itr+1) = var4_new;
        H = H_new;
        
        
    end
    
    
    if isempty(c_hat),
        c_hat = c_hat_new;
    else
        c_hat = c_hat_new*gamma_sync + c_hat*(1-gamma_sync);
    end
    
    
end



% Data decoding
llr_out = zeros(1, N_sym*N_dbps);

if (1),
    
    
    for m=1:N_sym,
        
        R = fft(rx_dec((1+m)*(N_gi+N_fft)+N_gi+1:(2+m)*(N_gi+N_fft)));
        R_eq = zeros(size(R));
        R_eq(P_active_mod) = R(P_active_mod) ./ H(P_active_mod);
        
        
        
        csi_ndc = 1;
        
        % slicer decisions
        if (enable_ndc),
            reiterate = 1;
            reiterate_count = 0;
            while reiterate,
                R_comp = R_eq;
                R_comp_best = R_eq;
                Err_best = 0;
                err = [];
                reiterate = 0;
                for it=1:length(alp),
                    
                    X_hat = zeros(size(R_comp));
                    
                    X_hat(P_data_mod) = slicer(R_comp(P_data_mod), N_bpscs);
                    X_hat(P_pilots_mod) = pilot_pattern * (2*pilot_polarity(m+1)-1);
                    
                    [d3, d5] = calcd(X_hat, N_fft);
                    U = [d3(P_active_mod) - T3*X_hat(P_active_mod); d5(P_active_mod) - T5*X_hat(P_active_mod)].';
                    A = diag(H(P_active_mod))' * diag(H(P_active_mod)); 
                    c_hat_new = inv(U'*A*U + lambda*eye(2))*U'*A*(R_eq(P_active_mod)-X_hat(P_active_mod)).';
                    
                    if isempty(c_hat),
                        c_hat = c_hat_new;
                    else
                        %err = [err; sum(abs(c_hat - c_hat_new).^2)];
                        %err = [err; sum(abs(R_eq(P_long)-X_hat(P_long)).^2)];
                        err = [err; c_hat_new(1)];
                        c_hat = c_hat_new*gamma_data + c_hat*(1-gamma_data);
                    end
                    
                    R_comp(P_active_mod) = R_eq(P_active_mod) - alp(it)*(U*c_hat).';
                    
                    if it==1,
                        err1 = sum(abs(R_eq(P_active_mod) - X_hat(P_active_mod)).^2);
                        Err_best = err1;
                    end
                    
                    err2 = sum(abs(R_comp(P_active_mod) - X_hat(P_active_mod)).^2);
                    if err2<Err_best,
                        R_comp_best(P_active_mod) = R_comp(P_active_mod);
                        Err_best = err2;
                    end
                    
                    if 0 && (it == length(alp)) && (mean(abs(err(1:end-1) - err(2:end))) > R_th),
                        c_hat = zeros(size(c_hat));
                        if reiterate_count<1,
                            reiterate = 1;
                            reiterate_count = reiterate_count + 1;
                            fprintf('R');
                        else
                            R_comp(P_active_mod) = R_eq(P_active_mod);
                            fprintf('X');
                        end
                    end
                    

                end
            end
            
            R_eq = R_comp_best;
            
          

            
        end
                
        data_rx = R_eq(P_data_mod);
        csi_rx = H(P_data_mod);
        llr = demapper_csi(data_rx, csi_rx*csi_ndc, N_bpscs);
        llr = depermuter_ac(llr, N_cbps, N_bpscs);
        llr_out((m-1)*N_cbps+1:m*N_cbps) = llr;
        
    end
end


llr_depunct = zeros(1, round(length(llr_out)*cr*2));

P_len = length(P_pattern);
P_nz = length(find(P_pattern == 1));
for k=1:length(llr_out)/P_nz,
    data_chunk = zeros(1,P_len);
    data_chunk(P_pattern==1) = llr_out(1+(k-1)*P_nz:k*P_nz);
    llr_depunct(1+(k-1)*P_len:k*P_len) = data_chunk;
end

% Remove pad bits before viterbi decoding
llr_depunct = llr_depunct(1:length(llr_depunct)-2*N_pad);

data_dec = vitdec(llr_depunct, trellis, 48, 'term', 'unquant');
% reinsert initial scrambler state for accurate BER at low SNR
data_dec(1:7) = [0 0 0 0 1 1 1]; 
data_descr = descrambler(data_dec);

% Remove tail bits after viterbi decoding
data_out = data_descr(N_SERVICE_FIELD+1:length(data_descr)-N_TAIL_FIELD);



end

