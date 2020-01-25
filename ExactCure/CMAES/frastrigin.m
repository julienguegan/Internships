function y = frastrigin(x)
    if size(x,1) < 2
        error('dimension must be greater one');
    end
    N = size(x,1);
    y = 10*N + sum(x(1:N).^2-10*cos(2*pi*x(1:N)));
end

