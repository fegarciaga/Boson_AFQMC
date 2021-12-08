function [phi_u, phi_d, O, w] = V_2(phi_u, phi_d, phi_T_u, phi_T_d, n_up, n_dn, O, w, Uab, N_par, deltau)
    %% In this function, the real Stratonovich-Hubbard transformation is implemented (NO phaseless condition)
    % aux1 will be used when correcting the overlap of permanents
    aux1=phi_T_u'*phi_u;
    v_u=sqrt(Uab*deltau/2)*(aux1*N_par/O-n_up);
    aux2=phi_T_d'*phi_d;
    v_d=sqrt(Uab*deltau/2)*(aux2*N_par/O-n_dn);
    %
    F=v_u-v_d;
    %% propagates the walker with the chosen auxiliary field
    aux_field=randn+F;
    phi_u=phi_u*exp(aux_field*sqrt(Uab/2*deltau));
    phi_d=phi_d*exp(-aux_field*sqrt(Uab/2*deltau));
    w=w*exp(-sqrt(Uab*deltau/2)*aux_field*(n_up-n_dn));
    % corrects overlap
    O_new=O-aux1-aux2+phi_T_u'*phi_u+phi_T_d'*phi_d;
    O_rat=O_new/O;
    x_path=max(0,cos(angle(O_rat^N_par)));
    w=real(w*exp(-sqrt(Uab*deltau/2)*aux_field*(n_up-n_dn))*O_rat^N_par*exp(F^2/2-aux_field*F))*x_path;
    O=O_new;
end