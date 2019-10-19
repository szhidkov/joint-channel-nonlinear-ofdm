function llr_out = demapper_csi(rx_in, H, mode)
%MAPPER Summary of this function goes here
%   Detailed explanation goes here

N0 = 1; % No scaling, regular Viterbi decoding
type = 1; % 0 - exact, 1 - approximate

if mode==1,
    % BPSK
    llr_out = -real(rx_in)/N0 .* abs(H).^2;
    
elseif mode==2,
    % QPSK
    llr_out = zeros(1, 2*length(rx_in));
    llr_out(1:2:end) = -real(rx_in)/N0 .* abs(H).^2;
    llr_out(2:2:end) = -imag(rx_in)/N0 .* abs(H).^2;
    
elseif mode==4,
    % 16-QAM
    q1 = 1 / sqrt(10);
    q3 = 3 / sqrt(10);
    
    llr_out = zeros(1,4*length(rx_in));
    
    for k=1:length(rx_in),
        re = real(rx_in(k));
        im = imag(rx_in(k));
        llr0 = min([(q3-re).^2 (q1-re).^2]) - min([(-q3-re).^2 (-q1-re).^2]) / N0;        
        llr1 = min([(q1-re).^2 (-q1-re).^2]) - min([(q3-re).^2 (-q3-re).^2]) / N0;
        llr2 = min([(q3-im).^2 (q1-im).^2]) - min([(-q3-im).^2 (-q1-im).^2]) / N0;        
        llr3 = min([(q1-im).^2 (-q1-im).^2]) - min([(q3-im).^2 (-q3-im).^2]) / N0;
        llr_out((k-1)*4+1) = llr0 * abs(H(k)).^2; 
        llr_out((k-1)*4+2) = llr1 * abs(H(k)).^2;  
        llr_out((k-1)*4+3) = llr2 * abs(H(k)).^2; 
        llr_out((k-1)*4+4) = llr3 * abs(H(k)).^2; 
    end
    
    
elseif mode==6,
    % 64-QAM
    
    q1 = 1 / sqrt(42);
    q3 = 3 / sqrt(42);
    q5 = 5 / sqrt(42);
    q7 = 7 / sqrt(42);
    
    llr_out = zeros(1,6*length(rx_in));
    
    for k=1:length(rx_in),
        re = real(rx_in(k));
        im = imag(rx_in(k));
        llr0 = min([(q7-re).^2 (q5-re).^2 (q3-re).^2 (q1-re).^2]) ...
             - min([(-q7-re).^2 (-q5-re).^2 (-q3-re).^2 (-q1-re).^2]) / N0;        
        llr1 = min([(q1-re).^2 (q3-re).^2 (-q3-re).^2 (-q1-re).^2]) ...
             - min([(-q7-re).^2 (-q5-re).^2 (q5-re).^2 (q7-re).^2]) / N0; 
        llr2 = min([(q5-re).^2 (q3-re).^2 (-q5-re).^2 (-q3-re).^2]) ...
             - min([(-q7-re).^2 (-q1-re).^2 (q7-re).^2 (q1-re).^2]) / N0;         
        llr3 = min([(q7-im).^2 (q5-im).^2 (q3-im).^2 (q1-im).^2]) ...
             - min([(-q7-im).^2 (-q5-im).^2 (-q3-im).^2 (-q1-im).^2]) / N0;        
        llr4 = min([(q1-im).^2 (q3-im).^2 (-q3-im).^2 (-q1-im).^2]) ...
             - min([(-q7-im).^2 (-q5-im).^2 (q5-im).^2 (q7-im).^2]) / N0; 
        llr5 = min([(q5-im).^2 (q3-im).^2 (-q5-im).^2 (-q3-im).^2]) ...
             - min([(-q7-im).^2 (-q1-im).^2 (q7-im).^2 (q1-im).^2]) / N0;  
         
        llr_out((k-1)*6+1) = llr0 * abs(H(k)).^2;  
        llr_out((k-1)*6+2) = llr1 * abs(H(k)).^2;  
        llr_out((k-1)*6+3) = llr2 * abs(H(k)).^2; 
        llr_out((k-1)*6+4) = llr3 * abs(H(k)).^2; 
        llr_out((k-1)*6+5) = llr4 * abs(H(k)).^2; 
        llr_out((k-1)*6+6) = llr5 * abs(H(k)).^2; 
        
    end   
 
elseif mode==8,
    % 256-QAM
    
    q1 = 1 / sqrt(170);
    q3 = 3 / sqrt(170);
    q5 = 5 / sqrt(170);
    q7 = 7 / sqrt(170);
    q9 = 9 / sqrt(170);
    q11 = 11 / sqrt(170);
    q13 = 13 / sqrt(170);
    q15 = 15 / sqrt(170);
    
    llr_out = zeros(1,8*length(rx_in));
    
    for k=1:length(rx_in),
        re = real(rx_in(k));
        im = imag(rx_in(k));
        
        llr0 = min([(q15-re).^2 (q13-re).^2 (q11-re).^2 (q9-re).^2 (q7-re).^2 (q5-re).^2 (q3-re).^2 (q1-re).^2]) ...
             - min([(-q15-re).^2 (-q13-re).^2 (-q11-re).^2 (-q9-re).^2 (-q7-re).^2 (-q5-re).^2 (-q3-re).^2 (-q1-re).^2]) / N0;     
         
        llr1 = min([(q7-re).^2 (q5-re).^2 (q3-re).^2 (q1-re).^2 (-q7-re).^2 (-q5-re).^2 (-q3-re).^2 (-q1-re).^2]) ...
             - min([(q15-re).^2 (q13-re).^2 (q11-re).^2 (q9-re).^2 (-q15-re).^2 (-q13-re).^2 (-q11-re).^2 (-q9-re).^2]) / N0; 
         
        llr2 = min([(q11-re).^2 (q9-re).^2 (q7-re).^2 (q5-re).^2 (-q11-re).^2 (-q9-re).^2 (-q7-re).^2 (-q5-re).^2]) ...
             - min([(q15-re).^2 (q13-re).^2 (q1-re).^2 (q3-re).^2 (-q15-re).^2 (-q13-re).^2 (-q1-re).^2 (-q3-re).^2]) / N0; 
   
        llr3 = min([(q13-re).^2 (q11-re).^2 (q5-re).^2 (q3-re).^2 (-q13-re).^2 (-q11-re).^2 (-q5-re).^2 (-q3-re).^2]) ...
             - min([(q15-re).^2 (q9-re).^2 (q7-re).^2 (q1-re).^2 (-q15-re).^2 (-q9-re).^2 (-q7-re).^2 (-q1-re).^2]) / N0; 

        llr4 = min([(q15-im).^2 (q13-im).^2 (q11-im).^2 (q9-im).^2 (q7-im).^2 (q5-im).^2 (q3-im).^2 (q1-im).^2]) ...
             - min([(-q15-im).^2 (-q13-im).^2 (-q11-im).^2 (-q9-im).^2 (-q7-im).^2 (-q5-im).^2 (-q3-im).^2 (-q1-im).^2]) / N0;     
         
        llr5 = min([(q7-im).^2 (q5-im).^2 (q3-im).^2 (q1-im).^2 (-q7-im).^2 (-q5-im).^2 (-q3-im).^2 (-q1-im).^2]) ...
             - min([(q15-im).^2 (q13-im).^2 (q11-im).^2 (q9-im).^2 (-q15-im).^2 (-q13-im).^2 (-q11-im).^2 (-q9-im).^2]) / N0; 
         
        llr6 = min([(q11-im).^2 (q9-im).^2 (q7-im).^2 (q5-im).^2 (-q11-im).^2 (-q9-im).^2 (-q7-im).^2 (-q5-im).^2]) ...
             - min([(q15-im).^2 (q13-im).^2 (q1-im).^2 (q3-im).^2 (-q15-im).^2 (-q13-im).^2 (-q1-im).^2 (-q3-im).^2]) / N0; 
   
        llr7 = min([(q13-im).^2 (q11-im).^2 (q5-im).^2 (q3-im).^2 (-q13-im).^2 (-q11-im).^2 (-q5-im).^2 (-q3-im).^2]) ...
             - min([(q15-im).^2 (q9-im).^2 (q7-im).^2 (q1-im).^2 (-q15-im).^2 (-q9-im).^2 (-q7-im).^2 (-q1-im).^2]) / N0; 
         
        llr_out((k-1)*8+1) = llr0 * abs(H(k)).^2;  
        llr_out((k-1)*8+2) = llr1 * abs(H(k)).^2;  
        llr_out((k-1)*8+3) = llr2 * abs(H(k)).^2; 
        llr_out((k-1)*8+4) = llr3 * abs(H(k)).^2; 
        llr_out((k-1)*8+5) = llr4 * abs(H(k)).^2; 
        llr_out((k-1)*8+6) = llr5 * abs(H(k)).^2; 
        llr_out((k-1)*8+7) = llr6 * abs(H(k)).^2; 
        llr_out((k-1)*8+8) = llr7 * abs(H(k)).^2; 
        
    end       
    
else
    % wrong QAM mode

end

