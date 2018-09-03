function [Rxx] = atc( x,timelag )
%finds Rxx(t) given x and t
N = length(x);
Rxx = sum((x(1:N-timelag)).*(x(timelag+1:N)));