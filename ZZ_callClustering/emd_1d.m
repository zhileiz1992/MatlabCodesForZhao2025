function d = emd_1d(p, q)
    d = sum(abs(cumsum(p) - cumsum(q)));
end
