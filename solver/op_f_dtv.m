function rhs = op_f_dtv (spv, msh, coeff)
  
 coeff = reshape (coeff, spv.ncomp, [], msh.nqn, msh.nel);

 rhs   = zeros (spv.ndof, 1);
 gradv  = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);
 ndir = size (gradv, 2);

 jacdet_weights = msh.jacdet .* msh.quad_weights;

 for iel = 1:msh.nel
   if (all (msh.jacdet(:,iel)))
     coeff_iel = reshape (coeff(:,:,:,iel), 1, msh.nqn, 1, 1);
     jacdet_iel = reshape (jacdet_weights(:, iel), [1, msh.nqn, 1, 1]);
     coeff_times_jw = bsxfun (@times, jacdet_iel, coeff_iel);

     gradv_iel = reshape (gradv(1,end,:,:,iel), 1, msh.nqn, spv.nsh_max, 1);

     aux_val = bsxfun (@times, coeff_times_jw, gradv_iel);
     rhs_loc = sum (sum (aux_val, 1), 2);

     indices = find (spv.connectivity(:,iel));
     rhs_loc = rhs_loc(indices); conn_iel = spv.connectivity(indices,iel);
     rhs(conn_iel) = rhs(conn_iel) + rhs_loc(:); 
   else
     warning ('geopdes:jacdet_zero_at_quad_node', 'op_f_gradv: singular map in element number %d', iel)
   end
 end
 
end
