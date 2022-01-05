function [E_ave,E_err, savedFileName]=PPMC_Bos(Lx,Ly,Lz,N_par,kx,ky,kz,U,Uab,jj,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_pc,itv_nrm,itv_Em,suffix)
tic; % start the  timer
bosonic_initialization; % initialize internal constants, form the trial wave function and assemble the initial population of walkers
format long;
flag_mea=0; %determine when a measurement should take place
E=0;
W=0;

% Preallocate arrays:
E_blk=zeros(N_blk,1); % array to store the energy measured in every block
W_blk=zeros(N_blk,1); % array to store the total weight in every block
%% Equilibration phase
p=parpool('local',2);
for i_blk=1:N_eqblk
    for j_step=1:N_blksteps
        [Phi, w, O, E, W] = bos_stepwlk(Phi, N_wlk, N_sites, w, O, n_up, n_dn, E, W, K_old, Proj_k_half, flag_mea, Phi_T, N_par, U, Uab, fac_norm, E_T, deltau);
        if mod(j_step,itv_pc)==0
            [Phi, w, O]=pop_cntrl_bos(Phi, w, O, N_wlk, N_sites); % population control
        end
        if mod(j_step,itv_nrm)==0
            [Phi, O, w]=bos_norm(Phi, Phi_T, N_wlk, O, w); % normalization control
        end
    end
end
%% Measurement phase    
for i_blk=1:N_blk
    for j_step=1:N_blksteps
        if mod(j_step,itv_Em)==0
            flag_mea=1;
        else
            flag_mea=0;
        end
        % propagate the walkers:
        [Phi, w, O, E_blk(i_blk), W_blk(i_blk)] = bos_stepwlk(Phi, N_wlk, N_sites, w, O, n_up, n_dn, E_blk(i_blk), W_blk(i_blk), K_old, Proj_k_half, flag_mea, Phi_T, N_par, U, Uab, fac_norm, E_T, deltau);
        if mod(j_step,itv_pc)==0
            [Phi, w, O]=pop_cntrl_bos(Phi, w, O, N_wlk, N_sites); % population control
        end
        if mod(j_step,itv_nrm)==0
            [Phi, O, w]=bos_norm(Phi, Phi_T, N_wlk, O, w); % normalization control
        end
        if mod(j_step, itv_Em)==0
            % update the exponent of the pre-factor exp(-deltau*(H-E_T))
            fac_norm=real(E_blk(i_blk)/W_blk(i_blk))*deltau+(-0.5*U*((n_up'*n_up)+(n_dn'*n_dn))-0.25*Uab*((n_up+n_dn)'*(n_up+n_dn))+0.25*Uab*((n_up-n_dn)'*(n_up-n_dn)))*deltau;
        end
    end
    E_blk(i_blk)=E_blk(i_blk)/W_blk(i_blk);
    display(strcat('E(',int2str(i_blk),')=',num2str(real(E_blk(i_blk)))))
end
delete(p);
%% Results
E=real(E_blk);
E_ave=mean(E)
E_err=std(E)/sqrt(N_blk)
% The total computational time:
time=toc() % stops the timer

%% Save data to a *.mat file
save (savedFileName, 'E', 'E_ave', 'E_err', 'time');
save (savedFileName, '-append', 'Lx', 'Ly','Lz', 'N_par', 'kx', 'ky','kz', 'U', 'tx', 'ty','tz');
save (savedFileName, '-append', 'deltau', 'N_wlk', 'N_blksteps', 'N_eqblk', 'N_blk', 'itv_pc','itv_nrm','itv_Em');
save (savedFileName, '-append', 'K', 'Phi_T');
