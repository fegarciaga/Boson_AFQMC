function [phi_u, phi_d, Aux_force, Log_W] = V_1(phi_u, phi_d, shift, n, Uab, deltau, Aux_force, Log_W)
    %% propagates the walker with the chosen auxiliary field
    exp_fact=randn+shift;
    Aux_force=Aux_force+shift*exp_fact;
    phi_u=phi_u*exp(exp_fact*sqrt(-Uab*deltau/2));
    phi_d=phi_d*exp(exp_fact*sqrt(-Uab*deltau/2));
    Log_W=Log_W+0.5*shift^2-shift*exp_fact+Uab*deltau/2*n*exp_fact;
end