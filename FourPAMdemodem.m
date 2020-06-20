function result = FourPAMdemodem(y,d)
result = zeros(length(y), 2);
result = char(result);
    for i = 1:(length(y))
        if(y(i) <= (-2*d))
            result(i,:) = "00";
        elseif((y(i) > (-2*d)) && (y(i) < 0))
            result(i,:) = "01";
        elseif((y(i) < 2*d) && (y(i) > 0))
            result(i,:) = "11";
        elseif(y(i) >= 2*d)
            result(i,:) = "10";
        end
    end  
end
