function [phi, w, O, E, W] = bos_stepwlk(phi, N_wlk, N_sites, w, O, n_up, n_dn, E, W, K, Proj_k_half, flag_mea, Phi_T, N_par, U, Uab, fac_norm, deltau)
    %% Propagate each walker:
    e=zeros(N_wlk,1); % Array containing the energy of each walker
    for i_wlk=1:N_wlk
        Phi=phi(:,i_wlk);
        if w(i_wlk)>0
             % multiply by the pre-factor exp(deltau*(E_T)) in the ground-state projector 
            w(i_wlk)=w(i_wlk)*exp(fac_norm);
            % propagate by the kinetic term exp(-1/2*deltau*K)
            [Phi, w(i_wlk), O(i_wlk)]=bos_halfK(Phi, w(i_wlk), O(i_wlk), Proj_k_half, Phi_T, N_par);
            if w(i_wlk)>0
                % propagate each lattice site of a walker by the potential term:
                for j_site=1:N_sites
                    if w(i_wlk)>0
                        [Phi(j_site), O(i_wlk), w(i_wlk)]=V_up_dn(Phi(j_site), Phi_T(j_site), n_up(j_site), O(i_wlk), w(i_wlk), U, N_par, deltau);
                    end
                end
%                 
                for j_site=1:N_sites
                    if w(i_wlk)>0
                       [Phi(j_site+N_sites), O(i_wlk), w(i_wlk)]=V_up_dn(Phi(j_site+N_sites), Phi_T(j_site+N_sites), n_dn(j_site), O(i_wlk), w(i_wlk), U, N_par, deltau);
                    end
                end
%                 
                for j_site=1:N_sites
                    if w(i_wlk)>0
                       [Phi(j_site), Phi(j_site+N_sites), O(i_wlk), w(i_wlk)]=V_1(Phi(j_site), Phi(j_site+N_sites), Phi_T(j_site), Phi_T(j_site+N_sites), n_up(j_site), n_dn(j_site), O(i_wlk), w(i_wlk), Uab, N_par, deltau);
                    end
                end
%                 
                for j_site=1:N_sites
                    if w(i_wlk)>0
                       [Phi(j_site), Phi(j_site+N_sites), O(i_wlk), w(i_wlk)]=V_2(Phi(j_site), Phi(j_site+N_sites), Phi_T(j_site), Phi_T(j_site+N_sites), n_up(j_site), n_dn(j_site), O(i_wlk), w(i_wlk), Uab, N_par, deltau);
                    end
                end
            end
            if w(i_wlk)>0
                % propagate by the kinetic term exp(-1/2*deltau*K)
                [Phi(:), w(i_wlk), O(i_wlk)]=bos_halfK(Phi(:), w(i_wlk), O(i_wlk), Proj_k_half, Phi_T, N_par);            
                if w(i_wlk)>0
                    % measure the energy if needed:
                    if flag_mea==1
                        [e(i_wlk)]=measure_b(K, Phi(:), Phi_T, O(i_wlk), N_sites, N_par, U, Uab);
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