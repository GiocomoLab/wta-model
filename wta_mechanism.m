function [output] = wta_mechanism(input,E)

% do Eqs 2 and 3 in the Lyttle paper
scaled_input = input - (1-E).*( ones(size(input,1),1)*max(input) );
scaled_input = scaled_input.*(scaled_input > 0);

% re-scale the firing rate maps, as per Lyttle paper
output = scaled_input./(ones(size(input,1),1)* max(scaled_input) );

return
