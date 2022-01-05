function [phi_u, phi_d, Log_W] = V_2(phi_u, phi_d, shift, n, Uab, deltau, Log_W)
    %% propagates the walker with the chosen auxiliary field
    exp_fact=randn+shift;
    phi_u=phi_u*exp(exp_fact*sqrt(Uab/2*deltau));
    phi_d=phi_d*exp(-exp_fact*sqrt(Uab/2*deltau));
    Log_W=Log_W+0.5*shift^2-shift*exp_fact-Uab*deltau/2*n*exp_fact;
end