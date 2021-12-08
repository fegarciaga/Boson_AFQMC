function [e] = measure_b(K, phi, Phi_T, O, N_sites, N_par, U, Uab)
    %% Calculate Green's function
    G=(Phi_T(:,1)*phi(:,1)')';
    M=diag(G');
    %%  Calculate energy due to one body operators:
    one_body_energy=sum(sum(K.'.*G))*N_par/O+U/2*sum(M)*N_par/O;
    %% Calculate energy due to two body operators:
    two_body_energy=(U*(M.'*M))*N_par*(N_par-1)/O^2;
    % Calculate the green's function taking into account the off diagonal terms
    L=M(1:N_sites);
    N=M(N_sites+1:2*N_sites);
    L_off=diag(G(1:N_sites,1+N_sites:2*N_sites)');
    N_off=diag(G(N_sites+1:2*N_sites,1:N_sites)');
    two_body_energy=two_body_energy+(Uab*(L.'*N+L_off.'*N_off))*N_par^2/O^2;
    %% calculate the total energy:
    e=one_body_energy+two_body_energy;
    if abs(O^2)<1e-3
%         display(phi);
%         display(abs(O));
        e=0;
    end
end