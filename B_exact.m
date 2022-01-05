function [e]= B_exact(U, Uab, t, jj)
    H=[U sqrt(2)*jj 0 -sqrt(2)*2*t 0 0 0 0 0 0; sqrt(2)*jj Uab sqrt(2)*jj 0 -2*t -2*t 0 0 0 0; 0 sqrt(2)*jj U 0 0 0 -sqrt(2)*2*t 0 0 0; -sqrt(2)*2*t 0 0 0 jj jj 0 -sqrt(2)*2*t 0 0; 0 -2*t 0 jj 0 0 jj 0 -2*t 0; 0 -2*t 0 jj 0 0 jj 0 -2*t 0; 0 0 -sqrt(2)*2*t 0 jj jj 0 0 0 -sqrt(2)*2*t; 0 0 0 -sqrt(2)*2*t 0 0 0 U sqrt(2)*jj 0; 0 0 0 0 -2*t -2*t 0 sqrt(2)*jj Uab sqrt(2)*jj; 0 0 0 0 0 0 -sqrt(2)*2*t 0 sqrt(2)*jj U];
    [psi,E] = eig(H);
    e=E(1);
end
