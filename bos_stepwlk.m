function [phi, w, O, E, W] = bos_stepwlk(phi, N_wlk, N_sites, w, O, n_up, n_dn, E, W, K, Proj_k_half, flag_mea, Phi_T, N_par, U, Uab, fac_norm, E_T, deltau)
    %% Propagate each walker:
    e=zeros(N_wlk,1); % Array containing the energy of each walker
    parfor i_wlk=1:N_wlk
        Phi=phi(:,i_wlk);
        if w(i_wlk)>0
            phase=N_par*angle(O(i_wlk));
            if cos(phase)<=0
                w(i_wlk)=0;
            end
            if w(i_wlk)>0
                n_up_t=N_par*diag((Phi_T(1:N_sites)').'*Phi(1:N_sites).')/O(i_wlk);
                n_dn_t=N_par*diag((Phi_T(1+N_sites:2*N_sites)').'*Phi(1+N_sites:2*N_sites).')/O(i_wlk);
                w(i_wlk)=w(i_wlk)*cos(phase);
                Aux_force=0;
                fact=measure_b(K, Phi(:), Phi_T, O(i_wlk), N_sites, N_par, U, Uab);
                % propagate by the kinetic term exp(-1/2*deltau*K)
                [Phi]=bos_halfK(Phi, Proj_k_half);
                % propagate each lattice site of a walker by the potential term:
                Log_W=-phase*sqrt(-1);
                for j_site=1:N_sites
                    [Phi(j_site), Aux_force, Log_W]=V_up_dn(Phi(j_site), sqrt(-1*U*deltau)*(n_up_t(j_site)-n_up(j_site)), n_up(j_site), U, deltau, Aux_force, Log_W);
                end
    %                 
                for j_site=1:N_sites
                    [Phi(j_site+N_sites), Aux_force, Log_W]=V_up_dn(Phi(j_site+N_sites), sqrt(-1*U*deltau)*(n_dn_t(j_site)-n_dn(j_site)), n_dn(j_site), U, deltau, Aux_force, Log_W);
                end
    %                 
                for j_site=1:N_sites
                    [Phi(j_site), Phi(j_site+N_sites), Aux_force, Log_W]=V_1(Phi(j_site), Phi(j_site+N_sites), sqrt(-1*Uab*deltau/2)*((n_up_t(j_site)+n_dn_t(j_site))-(n_up(j_site)+n_dn(j_site))), n_up(j_site)+n_dn(j_site), Uab, deltau, Aux_force, Log_W);
                end
    %                 
                for j_site=1:N_sites
                    [Phi(j_site), Phi(j_site+N_sites), Log_W]=V_2(Phi(j_site), Phi(j_site+N_sites), sqrt(Uab*deltau/2)*((n_up_t(j_site)-n_dn_t(j_site))-(n_up(j_site)-n_dn(j_site))), n_up(j_site)-n_dn(j_site), Uab, deltau, Log_W);
                end
                % propagate by the kinetic term exp(-1/2*deltau*K)
                [Phi]=bos_halfK(Phi, Proj_k_half); 
                % applies the phaseless condition 
                O_new=Phi_T'*Phi;
                fact2=measure_b(K, Phi(:), Phi_T, O_new, N_sites, N_par, U, Uab);
                w(i_wlk)=w(i_wlk)*exp(fac_norm);
                W_aux=exp(-deltau*0.5*(fact+fact2));
                x_path=max(0,cos(imag(Aux_force)));
                Log_W=Log_W+N_par*log(O_new/O(i_wlk));
                w(i_wlk)=w(i_wlk)*abs(exp(Log_W))*x_path;
                O(i_wlk)=O_new;
                if w(i_wlk)>0
                    % measure the energy if needed:
                    if flag_mea==1
                        [e(i_wlk)]=fact2;
                    end
                end   
            end
        end
        phi(:,i_wlk)=Phi;
    end
    %% Compute the ensemble's total energy and weight if measurement took place
    if flag_mea==1
        for i_wlk=1:N_wlk
            if w(i_wlk)>0
                E=E+e(i_wlk)*w(i_wlk);
                W=W+w(i_wlk);
            end
        end
    end
    
end