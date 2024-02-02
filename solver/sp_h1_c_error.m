function [errh1, errl2, errh1s, errh1_elem, errl2_elem, errh1s_elem] = sp_h1_c_error (sp, msh, u, uex, graduex, c)

grad_valu = sp_eval_msh (u, sp, msh, 'gradient');
grad_valu = reshape (grad_valu, sp.ncomp, msh.rdim, msh.nqn, msh.nel);

for idir = 1:msh.rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
end
grad_valex  = reshape (feval (graduex, x{:}), sp.ncomp, msh.rdim, msh.nqn, msh.nel);

c_val = reshape(c(x{:}),[1,1, msh.nqn, msh.nel]);

w = msh.quad_weights .* msh.jacdet;

[errl2, errl2_elem] = sp_l2_error (sp, msh, u, uex);
errh1s_elem = sum (reshape (sum (sum ((grad_valu - grad_valex).^2, 1).*c_val, 2), [msh.nqn, msh.nel]) .* w);
errh1s = sqrt (sum (errh1s_elem));

errh1  = sqrt (errl2^2 + errh1s^2);

errh1_elem  = sqrt (errl2_elem.^2 + errh1s_elem);
errh1s_elem = sqrt (errh1s_elem);

end
