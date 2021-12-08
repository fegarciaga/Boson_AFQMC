function [phi, O, w] = V_up_dn(phi, phi_T, n, O, w, U, N_par, deltau)
    %% Calculates optimal shift in the auxiliary field sample
    % aux will be used when correcting the overlap of permanents
    aux=phi_T'*phi;
    F=-sqrt(-1)*sqrt(deltau*U)*(aux*N_par/O-n);
    %% propagates the walker with the chosen auxiliary field
    aux_field=randn+F;
    phi=phi*exp(-sqrt(-1)*aux_field*sqrt(U*deltau));
    w=w*exp(sqrt(-1)*aux_field*sqrt(U*deltau)*n);
    % corrects overlap
    O_new=O-aux+phi_T'*phi;
    O_rat=O_new/O;
    x_path=max(0,cos(angle(exp(sqrt(-1)*aux_field*sqrt(U*deltau)*n)*(O_rat^N_par))));
    w=real(w*O_rat^N_par*exp(F^2/2-aux_field*F))*x_path;
    O=O_new;
end