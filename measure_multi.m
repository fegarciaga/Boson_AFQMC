function [e]=measure_multi(K, phi, w, N_sites, N_par, N_per, U, Uab)
    e=0;
    for i=1:N_per
        for j=1:N_per
            O=phi(:,i)'*phi(:,j);
            e_loc=measure_b(K, phi(:,i), phi(:,j), O, N_sites, N_par, U, Uab);
            e=e+w(i)*w(j)'*e_loc;
        end
    end
end