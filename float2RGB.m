function [ r,g,b ] = float2RGB( value )
% 将0-1之间的数值转换成RGB
    if value < 0
        value = 0;
    elseif value > 1
        value = 1;
    end
    
    if value < 1/6     % 000-001
        r = 0;
        g = 0;
        b = value*6;
    elseif value < 1/3 % 001-011
        r = 0;
        g = 6*value-1;
        b = 1;
    elseif value < 1/2 % 011-010
        r = 0;
        g = 1;
        b = 3-6*value;
    elseif value < 2/3 % 010-110
        r = 6*value-3;
        g = 1;
        b = 0;
    elseif value < 5/6 % 110-100
        r = 1;
        g = 5-6*value;
        b = 0;
    else               % 100-101
        r = 1;
        g = 0;
        b = 6*value-5;
    end
end

