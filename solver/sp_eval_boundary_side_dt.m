function sp_side = sp_eval_boundary_side_dt (sp, msh_side)

  iside = msh_side.side_number;
  sp_side = sp_precompute (sp.boundary(iside), msh_side, 'value', false,'gradient',true);

end