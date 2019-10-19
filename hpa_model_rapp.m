function [ tx_out ] = hpa_model_rapp( tx_in, A_sat, v )
%HPA_MODEL Summary of this function goes here
%   Detailed explanation goes here

ro = abs(tx_in);
phi = exp(1i*angle(tx_in)); %tx_in ./ (ro+0.0000000001);
tx_out = ro.*(1+(ro/A_sat).^v).^(-1/v) .* phi;


end

