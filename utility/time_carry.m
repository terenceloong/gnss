function t = time_carry(t)
% Î¢Ãë¡¢ºÁÃë¡¢Ãë½øÎ»

if t(3)>=1000
    t(2) = t(2) + 1;
    t(3) = t(3) - 1000;
elseif t(3)<0
    t(2) = t(2) - 1;
    t(3) = t(3) + 1000;
end
if t(2)>=1000
    t(1) = t(1) + 1;
    t(2) = t(2) - 1000;
elseif t(2)<0
    t(1) = t(1) - 1;
    t(2) = t(2) + 1000;
end

end