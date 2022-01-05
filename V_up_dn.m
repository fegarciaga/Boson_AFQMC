function [phi, Aux_force, Log_W] = V_up_dn(phi, shift, n, U, deltau, Aux_force, Log_W)
    %% propagates the walker with the chosen auxiliary field
    exp_fact=randn+shift;
    Aux_force=Aux_force+shift*exp_fact;
    phi=phi*exp(exp_fact*sqrt(-U*deltau));
    Log_W=Log_W+0.5*shift^2-shift*exp_fact+U*deltau*n*exp_fact;
end