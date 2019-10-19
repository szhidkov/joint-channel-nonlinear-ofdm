function data_out = permuter_ac(data_in, Ncbps, Nbpsc)
%PERMUTER Summary of this function goes here
%   Detailed explanation goes here

s = max([Nbpsc/2 1]);
k = 0:1:Ncbps/2-1;

% Segment parser
Nes = 1;
i0 = 2*s*Nes*floor(k/(s*Nes)) + 0*s*Nes + mod(k, s*Nes) + 1;
i1 = 2*s*Nes*floor(k/(s*Nes)) + 1*s*Nes + mod(k, s*Nes) + 1;

data_in_0 = data_in(i0);
data_in_1 = data_in(i1);

% BCC interleaver (two interleavers) 
N_col = 26;
N_row = Nbpsc*9;

perm1 = N_row*mod(k, N_col) + floor(k/N_col) + 1;
perm2 = s*floor(k/s) + mod(k+Ncbps/2-floor(N_col*k/(Ncbps/2)), s) + 1;

data_out_0 = zeros(size(data_in_0));
data_out_1 = zeros(size(data_in_1));

data_out_0(perm1) = data_in_0;
data_out_0(perm2) = data_out_0;

data_out_1(perm1) = data_in_1;
data_out_1(perm2) = data_out_1;

% Segment deparser
data_out = zeros(size(data_in));
data_out(i0) = data_out_0;
data_out(i1) = data_out_1;


end

