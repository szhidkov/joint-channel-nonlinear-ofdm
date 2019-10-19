function tx_out = mapper(data, mode)
%MAPPER Summary of this function goes here
%   Detailed explanation goes here
if mode==1,
    % BPSK
    tx_out = 2*data - 1;
    
elseif mode==2,
    % QPSK
    tx_out = (2*data(1:2:end)-1) + 1i*(2*data(2:2:end)-1);
    tx_out = tx_out*(1/sqrt(2));
    
elseif mode==4,
    % 16-QAM
    b0 = data(1:4:end);
    b1 = data(2:4:end);
    b2 = data(3:4:end);
    b3 = data(4:4:end);    
    tx_out = zeros(1,length(data)/4);
    for k=1:length(tx_out),
       if b0(k)==0 && b1(k)==0,
           tx_out(k) = -3;
       elseif b0(k)==0 && b1(k)==1,
           tx_out(k) = -1;
       elseif b0(k)==1 && b1(k)==1,
           tx_out(k) = 1;
       elseif b0(k)==1 && b1(k)==0,
           tx_out(k) = 3;
       end
       if b2(k)==0 && b3(k)==0,
           tx_out(k) = tx_out(k) - 1i*(3);
       elseif b2(k)==0 && b3(k)==1,
           tx_out(k) = tx_out(k) - 1i*(1);
       elseif b2(k)==1 && b3(k)==1,
           tx_out(k) = tx_out(k) + 1i*(1);
       elseif b2(k)==1 && b3(k)==0,
           tx_out(k) = tx_out(k) + 1i*(3);
       end       
    end
    tx_out = tx_out * (1/sqrt(10));
    
elseif mode==6,
    % 64-QAM
    b0 = data(1:6:end);
    b1 = data(2:6:end);
    b2 = data(3:6:end);
    b3 = data(4:6:end);
    b4 = data(5:6:end);
    b5 = data(6:6:end);
    
    tx_out = zeros(1,length(data)/6);
    
    for k=1:length(tx_out),
        
       if b0(k)==0 && b1(k)==0 && b2(k)==0,
           tx_out(k) = -7;
       elseif b0(k)==0 && b1(k)==0 && b2(k)==1,
           tx_out(k) = -5;
       elseif b0(k)==0 && b1(k)==1 && b2(k)==1,
           tx_out(k) = -3;
       elseif b0(k)==0 && b1(k)==1 && b2(k)==0,
           tx_out(k) = -1;
       elseif b0(k)==1 && b1(k)==1 && b2(k)==0,
           tx_out(k) = 1;           
       elseif b0(k)==1 && b1(k)==1 && b2(k)==1,
           tx_out(k) = 3;
       elseif b0(k)==1 && b1(k)==0 && b2(k)==1,
           tx_out(k) = 5;        
       elseif b0(k)==1 && b1(k)==0 && b2(k)==0,
           tx_out(k) = 7;              
       end
       
       if b3(k)==0 && b4(k)==0 && b5(k)==0,
           tx_out(k) = tx_out(k) - 1i*(7);
       elseif b3(k)==0 && b4(k)==0 && b5(k)==1,
           tx_out(k) = tx_out(k) - 1i*(5);
       elseif b3(k)==0 && b4(k)==1 && b5(k)==1,
           tx_out(k) = tx_out(k) - 1i*(3);
       elseif b3(k)==0 && b4(k)==1 && b5(k)==0,
           tx_out(k) = tx_out(k) - 1i*(1);
       elseif b3(k)==1 && b4(k)==1 && b5(k)==0,
           tx_out(k) = tx_out(k) + 1i*(1);        
       elseif b3(k)==1 && b4(k)==1 && b5(k)==1,
           tx_out(k) = tx_out(k) + 1i*(3);  
       elseif b3(k)==1 && b4(k)==0 && b5(k)==1,
           tx_out(k) = tx_out(k) + 1i*(5);        
       elseif b3(k)==1 && b4(k)==0 && b5(k)==0,
           tx_out(k) = tx_out(k) + 1i*(7);             
       end
       
    end
    tx_out = tx_out * (1/sqrt(42));    
    
elseif mode==8,
    % 256-QAM
    b0 = data(1:8:end);
    b1 = data(2:8:end);
    b2 = data(3:8:end);
    b3 = data(4:8:end);
    b4 = data(5:8:end);
    b5 = data(6:8:end);
    b6 = data(7:8:end);
    b7 = data(8:8:end);
    
    tx_out = zeros(1,length(data)/8);   

    for k=1:length(tx_out),
        
       if b0(k)==0 && b1(k)==0 && b2(k)==0 && b3(k)==0,
           tx_out(k) = -15;
       elseif b0(k)==0 && b1(k)==0 && b2(k)==0 && b3(k)==1,
           tx_out(k) = -13;
       elseif b0(k)==0 && b1(k)==0 && b2(k)==1 && b3(k)==1,
           tx_out(k) = -11;
       elseif b0(k)==0 && b1(k)==0 && b2(k)==1 && b3(k)==0,
           tx_out(k) = -9;
       elseif b0(k)==0 && b1(k)==1 && b2(k)==1 && b3(k)==0,
           tx_out(k) = -7;
       elseif b0(k)==0 && b1(k)==1 && b2(k)==1 && b3(k)==1,
           tx_out(k) = -5;    
       elseif b0(k)==0 && b1(k)==1 && b2(k)==0 && b3(k)==1,
           tx_out(k) = -3;
       elseif b0(k)==0 && b1(k)==1 && b2(k)==0 && b3(k)==0,
           tx_out(k) = -1;
       elseif b0(k)==1 && b1(k)==0 && b2(k)==0 && b3(k)==0,
           tx_out(k) = 15;
       elseif b0(k)==1 && b1(k)==0 && b2(k)==0 && b3(k)==1,
           tx_out(k) = 13;
       elseif b0(k)==1 && b1(k)==0 && b2(k)==1 && b3(k)==1,
           tx_out(k) = 11;
       elseif b0(k)==1 && b1(k)==0 && b2(k)==1 && b3(k)==0,
           tx_out(k) = 9;
       elseif b0(k)==1 && b1(k)==1 && b2(k)==1 && b3(k)==0,
           tx_out(k) = 7;
       elseif b0(k)==1 && b1(k)==1 && b2(k)==1 && b3(k)==1,
           tx_out(k) = 5;    
       elseif b0(k)==1 && b1(k)==1 && b2(k)==0 && b3(k)==1,
           tx_out(k) = 3;
       elseif b0(k)==1 && b1(k)==1 && b2(k)==0 && b3(k)==0,
           tx_out(k) = 1;           
       end
       
       if b4(k)==0 && b5(k)==0 && b6(k)==0 && b7(k)==0,
           tx_out(k) = tx_out(k) + 1i*(-15);
       elseif b4(k)==0 && b5(k)==0 && b6(k)==0 && b7(k)==1,
           tx_out(k) = tx_out(k) + 1i*(-13);
       elseif b4(k)==0 && b5(k)==0 && b6(k)==1 && b7(k)==1,
           tx_out(k) = tx_out(k) + 1i*(-11);
       elseif b4(k)==0 && b5(k)==0 && b6(k)==1 && b7(k)==0,
           tx_out(k) = tx_out(k) + 1i*(-9);
       elseif b4(k)==0 && b5(k)==1 && b6(k)==1 && b7(k)==0,
           tx_out(k) = tx_out(k) + 1i*(-7);
       elseif b4(k)==0 && b5(k)==1 && b6(k)==1 && b7(k)==1,
           tx_out(k) = tx_out(k) + 1i*(-5);    
       elseif b4(k)==0 && b5(k)==1 && b6(k)==0 && b7(k)==1,
           tx_out(k) = tx_out(k) + 1i*(-3);
       elseif b4(k)==0 && b5(k)==1 && b6(k)==0 && b7(k)==0,
           tx_out(k) = tx_out(k) + 1i*(-1);
       elseif b4(k)==1 && b5(k)==0 && b6(k)==0 && b7(k)==0,
           tx_out(k) = tx_out(k) + 1i*(15);
       elseif b4(k)==1 && b5(k)==0 && b6(k)==0 && b7(k)==1,
           tx_out(k) = tx_out(k) + 1i*(13);
       elseif b4(k)==1 && b5(k)==0 && b6(k)==1 && b7(k)==1,
           tx_out(k) = tx_out(k) + 1i*(11);
       elseif b4(k)==1 && b5(k)==0 && b6(k)==1 && b7(k)==0,
           tx_out(k) = tx_out(k) + 1i*(9);
       elseif b4(k)==1 && b5(k)==1 && b6(k)==1 && b7(k)==0,
           tx_out(k) = tx_out(k) + 1i*(7);
       elseif b4(k)==1 && b5(k)==1 && b6(k)==1 && b7(k)==1,
           tx_out(k) = tx_out(k) + 1i*(5);    
       elseif b4(k)==1 && b5(k)==1 && b6(k)==0 && b7(k)==1,
           tx_out(k) = tx_out(k) + 1i*(3);
       elseif b4(k)==1 && b5(k)==1 && b6(k)==0 && b7(k)==0,
           tx_out(k) = tx_out(k) + 1i*(1);          
       end
       
    end
    tx_out = tx_out * (1/sqrt(170));       
    
else
    % wrong QAM mode

end

