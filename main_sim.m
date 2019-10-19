%
%  Simulation of joint nonlinear distortion compensation 
%  and channel estimation in OFDM receiver (Bussgang-type) 
%  simulation parameters are set to replicate WiFi system (ac)
%  
%  For more details on the algorithm see the paper:
%  S. V. Zhidkov, "Joint Channel Estimation and Nonlinear 
%  Distortion Compensation in OFDM Receivers," Arxiv, 2016
%  https://arxiv.org/pdf/1612.09222v1 
%  and presentation at: https://cifrasoft.com/people/szhidkov/
%  
%  Copyright (c) 2016-2017 Sergey Zhidkov
%  Licensed under the MIT license.
%  See LICENSE file in the project root for details. 
%

% Main simulation parameters
PSDU_LENGTH = 400;       % Length of PSDU (symbols) in TX packet
MSC_index = 8;           % Modulation and Coding Scheme index 
N_packets = 100;         % Number of TX packets to simulate per trial
channel_model = 'AWGN';  % Channel model: 'AWGN', 'A', 'B', 'C'
A_sat = 2.7;             % PA saturation level, see: hpa_model_rapp()

% Array of SNR values for simulation 
SNR = [17.5 18.75 20 21.25 22.5 23.75 25 26.25 27.5 30 35]; % for MSC=8
%SNR = [12.5 15 17.5 18.75 20 21.25 22.5 23.75 25 26.25 30]; % for MSC=6 

% Overriden parameters
enable_hpa = 1;          % Enable power amplifier switch (overriden)
enable_ndc = 1;          % Enable nonlinear distortion compensation (overriden)
enable_fec = 1;          % Enable FEC feedback (overriden)

% System parameters
Fs = 160;  % Only 160Mhz mode is currently supported
L = 5;     % Oversampling factor

min_packet_errors = 20; 
min_packets = 200;

BER1 = zeros(size(SNR));
PER1 = zeros(size(SNR));
BER2 = zeros(size(SNR));
PER2 = zeros(size(SNR));
BER3 = zeros(size(SNR));
PER3 = zeros(size(SNR));
BER_lin = zeros(size(SNR));
PER_lin = zeros(size(SNR));
BER_nln = zeros(size(SNR));
PER_nln = zeros(size(SNR));

for u=1:length(SNR)
    
    % 5 trials per SNR value:
    % 1: HPA - on, NDC - on, 1 iteration
    % 2: HPA - on, NDC - on, 2 iterations
    % 3: HPA - on, NDC - on, 3 iterations
    % 4: HPA - on, NDC - off (conventional receiver)
    % 5: HPA - off, NDC - off (linear channel)
    for trial_no=1:5
        
        if trial_no==1
            enable_hpa = 1;
            enable_ndc = 1;
            enable_fec = 0;
            n_iter = 1;
        elseif trial_no==2
            enable_hpa = 1;
            enable_ndc = 1;
            enable_fec = 0;
            n_iter = 2;
        elseif trial_no==3
            enable_hpa = 1;
            enable_ndc = 1;
            enable_fec = 0;
            n_iter = 3;
        elseif trial_no==4
            enable_hpa = 0;
            enable_ndc = 0;
            enable_fec = 0;
            n_iter = 1;
        elseif trial_no==5
            enable_hpa = 1;
            enable_ndc = 0;
            enable_fec = 0;
            n_iter = 1;    
        end
        
        c_hat = [];
        
        BER = 0;
        PER = 0;
        
        N_packets_final = N_packets;
        
        obo = 0;

        for pak = 1:N_packets
            
            data_psdu = randi([0 1], 1, PSDU_LENGTH*8);
            
            tx_out = tx_burst_ac(data_psdu, MSC_index, 160);
            
            tx_upsample = resample(tx_out, L, 1, 32); 
            

            if enable_hpa == 1
                A_sat_norm = A_sat*std(tx_upsample);
                tx_nln = hpa_model_rapp(tx_upsample, A_sat_norm, 2);
                
                %tx_nln = tx_upsample;
                
                obo = obo + 10*log10(A_sat_norm.^2 / mean(abs(tx_nln).^2));
       
                % Uncomment line below to plot spectrum with specral mask
                %f_msk = [0 79 81 160 240 Fs*L/2 Fs*L-240 Fs*L-160 Fs*L-81 Fs*L-79 Fs*L-0];
                %p_msk = [0 0 -20 -28 -40 -40 -40 -28 -20 0 0];
                %plot(f_msk, p_msk, 'r');
                %hold on
                %psd(tx_nln/std(tx_nln)/sqrt(6.5), 1024, Fs*L)
                %stop
                
            else
                tx_nln = tx_upsample;
            end
            

         
            % Multipath channel
            tx_ch = multipath_channel_ac(tx_nln, channel_model, SNR(u), L);
            %tx_ch = tx_nln;
            
                  
            % RX front end
            rx_in = resample(tx_ch, 2, L, 10);
            c_hat = [];
            [data_out, c_hat] = rx_demod_ac(rx_in, MSC_index, Fs, PSDU_LENGTH, enable_ndc, c_hat, n_iter);

            
            ber_pak = sum((data_psdu - data_out)~=0) / length(data_psdu);
            per_pak = (ber_pak>0)*1;
            
                  
            
            if per_pak>0
                fprintf('*');
            else
                fprintf('.');
            end
            
            
            
            BER = BER + ber_pak;
            PER = PER + per_pak;
            
            % Check if we already get half of expected packet errors
            % and let our simulation run twice as long as till now
            if (pak>=min_packets/2) && (PER > min_packet_errors/2) && (N_packets_final==N_packets), 
                N_packets_final = 2*pak;
            end
            
            if mod(pak,60)==0
                fprintf('\n');
            end
            
            if pak >= N_packets_final
                break;
            end
            
        end
        
        fprintf('\n');
        
        
        N_packets_final = min([N_packets N_packets_final]);
        
        BER = BER / N_packets_final;
        PER = PER / N_packets_final;
        
        fprintf('Trial: %i | SNR=%2.2f | BER=%2.8f | PER=%2.8f | OBO = %2.8f\n', trial_no, SNR(u), BER, PER, obo/N_packets_final);
        
        if trial_no==1
            BER1(u) = BER;
            PER1(u) = PER;
        elseif trial_no==2
            BER2(u) = BER;
            PER2(u) = PER;
        elseif trial_no==3
            BER3(u) = BER;
            PER3(u) = PER; 
        elseif trial_no==4
            BER_lin(u) = BER;
            PER_lin(u) = PER;
        elseif trial_no==5
            BER_nln(u) = BER;
            PER_nln(u) = PER;            
        end
        
    end
    
    semilogy(SNR, BER_lin, 'b-o');
    hold on
    semilogy(SNR, PER_lin, 'r--o');
    semilogy(SNR, BER_nln, 'b-v');
    semilogy(SNR, PER_nln, 'r--v');
    semilogy(SNR, BER1, 'b-x');
    semilogy(SNR, PER1, 'r--x');
    semilogy(SNR, BER2, 'b-x');
    semilogy(SNR, PER2, 'r--x');
    semilogy(SNR, BER3, 'b-x');
    semilogy(SNR, PER3, 'r--x'); 
    hold off
    pause(0.1)
    
end
 
% Save simulation results in mat-file
save results.mat SNR BER_lin PER_lin BER_nln PER_nln BER1 PER1 BER2 PER2 BER3 PER3 A_sat N_packets MSC_index channel_model

