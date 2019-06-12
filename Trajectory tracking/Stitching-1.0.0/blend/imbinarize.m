function out = imbinarize(in)
    sumin = sum(in, 3);
    out = sumin;
    out(sumin == 0) = 0;
    out(sumin ~= 0) = 1;
end

