function uHil = transform(u)

N = size(u,1);
M = size(u,2);
uHil = zeros(N,M);
for i = 1:N
    sigt = imag(hilbert(u(i,:)));
    uHil(i,:) = sigt;
end


end