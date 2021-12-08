function [phi, w, O] = bos_halfK(phi, w, O, Proj_k_half, Phi_T, N_par)
    %% propagate the walker by exp(-deltau*K/2)
    phi=Proj_k_half*phi;
    %% update the overlap
    O_new=Phi_T(:)'*phi(:);
    % calculate the new overlap
    O_ratio=O_new/O;
    %% enforce the constrained path condition
    % If the new weight is negative (O_ratio<0), kill the walker by setting its weight to zero
    if(O_ratio<0)
        w=0;
    else
        % real(O_ratio) enforces the phase-free approximation in case of complex phases (because the condition O_ratio>0 only checks the real part of O_ratio)
        O=O_new;
        x_path=max(0,cos(angle(O_ratio^N_par)));
        w=real(w*O_ratio^N_par)*x_path;
    end
end