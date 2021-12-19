%% 
function value = circ(i,max)
% Circ wraps any integer index around a max value.  Creates quick circular
% indices
%   Detailed explanation goes here
value = i;
while( (value > max) || (value < 1) )
    if(value > max)
        value = value - max;
    elseif(value < 1)
        value = value + max;
    end
end

end

