function [ a ] = levinson_recursion( x, N )
% Finds aN given x and N
Rxx = zeros(1,N+1);
for i = 1:N+1
    Rxx(i) = atc(x,i-1);
end
a = zeros(N,1);

a(1,1) = (-1*(atc(x,1)))/atc(x,0);
for i = 1 : N-1
    error = Rxx(1) + sum(a(1:i).*transpose((Rxx(2:i+1))));
    alpha = Rxx(i+2) + sum(a(1:i).*transpose((fliplr(Rxx(2:i+1)))));
    k = (-1*alpha)/error;
    a(1:i) = a(1:i) + k*(flipud(a(1:i)));
    a(i+1) = k;
end

