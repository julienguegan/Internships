% Log10 of tail of binomial distribution by Hoeffding approximation

function y = log_bin(n, k, p)

r =  k / n;
if( r <= p )
    y = 0;
elseif( n == k)
    y = k * log10(p);
else
    y = k * log10(p/r) + (n-k) * log10( (1-p)/(1-r) );
end

end