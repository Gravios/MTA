function cfr(range)
if (numel(range) == 1)  &&  (range > 0)
    for n = 1:range
        close(figure(n));
    end
else
    for n = range(:)'
        close(figure(n));
    end
end

        