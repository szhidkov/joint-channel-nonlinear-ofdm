function rx_out = slicer(rx_in, mode)
%MAPPER Summary of this function goes here
%   Detailed explanation goes here

N0 = 1; % No scaling, regular Viterbi decoding
type = 1; % 0 - exact, 1 - approximate

if mode==1,
    
    % BPSK
    rx_out = sign(real(rx_in)) * 1;
    
elseif mode==2,
    
    % QPSK
    rx_out = (sign(real(rx_in)) + 1i*sign(imag(rx_in))) / sqrt(2);
    
elseif mode==4,
    
    % 16-QAM
    rx_out = zeros(size(rx_in)); 
        
    for k=1:length(rx_in),
        
        re = real(rx_in(k)) * sqrt(10);
        im = imag(rx_in(k)) * sqrt(10);
        
        if re<-2,
            re = -3;
        elseif re<0,
            re = -1;
        elseif re<2,
            re = 1;
        else
            re = 3;
        end
 
        if im<-2,
            im = -3;
        elseif im<0,
            im = -1;
        elseif im<2,
            im = 1;
        else
            im = 3;
        end
        
        rx_out(k) = (re + 1i*im) / sqrt(10);
        
    end   

    
    
elseif mode==6,
    
    % 64-QAM
    rx_out = zeros(size(rx_in)); 
        
    for k=1:length(rx_in),
        
        re = real(rx_in(k)) * sqrt(42);
        im = imag(rx_in(k)) * sqrt(42);
        
        if re<-6,
            re = -7;
        elseif re<-4,
            re = -5;
        elseif re<-2,
            re = -3;
        elseif re<0
            re = -1;
        elseif re<2,
            re = 1;
        elseif re<4,
            re = 3;
        elseif re<6
            re = 5;
        else 
            re = 7; 
        end

        if im<-6,
            im = -7;
        elseif im<-4,
            im = -5;
        elseif im<-2,
            im = -3;
        elseif im<0
            im = -1;
        elseif im<2,
            im = 1;
        elseif im<4,
            im = 3;
        elseif im<6
            im = 5;
        else 
            im = 7; 
        end

        rx_out(k) = (re + 1i*im) / sqrt(42);
        
    end   
	
elseif mode==8,
    
    % 256-QAM
    rx_out = zeros(size(rx_in)); 
        
    for k=1:length(rx_in),
        
        re = real(rx_in(k)) * sqrt(170);
        im = imag(rx_in(k)) * sqrt(170);

        if re<-14,
            re = -15;
        elseif re<-12,
            re = -13;
        elseif re<-10,
            re = -11;
        elseif re<-8,
            re = -9;
        elseif re<-6,
            re = -7;
        elseif re<-4,
            re = -5;
        elseif re<-2,
            re = -3;        
        elseif re<0,
            re = -1;
        elseif re<2,
            re = 1;
        elseif re<4,
            re = 3;
        elseif re<6,
            re = 5;
        elseif re<8,
            re = 7;
        elseif re<10,
            re = 9;
        elseif re<12,
            re = 11;
        elseif re<14, 
            re = 13;
		else
			re = 15;
        end

        if im<-14,
            im = -15;
        elseif im<-12,
            im = -13;
        elseif im<-10,
            im = -11;
        elseif im<-8,
            im = -9;
        elseif im<-6,
            im = -7;
        elseif im<-4,
            im = -5;
        elseif im<-2,
            im = -3;        
        elseif im<0,
            im = -1;
        elseif im<2,
            im = 1;
        elseif im<4,
            im = 3;
        elseif im<6,
            im = 5;
        elseif im<8,
            im = 7;
        elseif im<10,
            im = 9;
        elseif im<12,
            im = 11;
        elseif im<14, 
            im = 13;
		else
			im = 15;
        end

        rx_out(k) = (re + 1i*im) / sqrt(170);
        
    end   
    	
    
else
    % wrong QAM mode
	error('Wrong QAM-mode');
end

