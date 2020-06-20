function result = FourPAMmodem(x,d)%innput is array
    result = zeros(length(x),1);
    for i = 1:(length(x))
        if(x(i,:) == '00')
            result(i) = -3*d;%(2*1-1-4)*1; %-3d
        elseif(x(i,:)== '01')
            result(i) = -1*d;%(2*2-1-4)*1; %-1d
        elseif(x(i,:)== '11')
            result(i) = 1*d;%(2*3-1-4)*1; %1d
        else %x = "11"
            result(i) = 3*d;%(2*4-1-4)*1; %3d
        end
    end
end