function scrambled_data = descrambler(data)
%DESCRAMBLER Summary of this function goes here
%   Detailed explanation goes here

initial_state = data(7:-1:1);

state = initial_state;
scrambled_data = zeros(size(data));

for k=8:length(data),
   s_bit = mod(state(7) + state(4), 2);
   scrambled_data(k) = mod(data(k) + s_bit, 2);
   state = [s_bit state(1:6)];
end




