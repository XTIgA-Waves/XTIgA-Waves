function varargout = op_mat_hostab_tp (space1, space2, msh,shape_functions_der1,shape_functions_der2,der)

A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);

    sp1_col = ho_sp_evaluate_col (space1, msh_col,shape_functions_der1);
    %sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, 'gradient', true, 'laplacian', false);
    
    sp2_col = ho_sp_evaluate_col (space2, msh_col,shape_functions_der2);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    A = A + op_mat_hostab (sp1_col, sp2_col, msh_col,der);
  end
   
  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
