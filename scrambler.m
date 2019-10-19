function scrambled_data = scrambler(data, initial_state)
%SCRAMBLER Summary of this function goes here
%   Detailed explanation goes here

if length(initial_state) ~= 7,
   error('Initial state should consists of 7 bit') 
end

state = initial_state;
scrambled_data = zeros(size(data));

for k=1:length(data),
   s_bit = mod(state(7) + state(4), 2);
   scrambled_data(k) = mod(data(k) + s_bit, 2);
   state = [s_bit state(1:6)];
end




