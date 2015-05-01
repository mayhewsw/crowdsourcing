function [ pHat, tHat ] = mv( A )
pHat = 0;
tHat = sign(sum(A, 2));