function sp = ho_sp_evaluate_col_param (space,msh,shape_functions_der)

sp_univ = space.sp_univ;

elem_list{1} = msh.colnum;
for idim = 2:msh.ndim
  elem_list{idim} = 1:msh.nel_dir(idim);
end

for idim = 1:msh.ndim
  nsh_dim{idim} = sp_univ(idim).nsh(elem_list{idim});
end

[nsh_grid{1:msh.ndim}] = ndgrid (nsh_dim{:});
nsh = 1;
for idim = 1:msh.ndim
  nsh = nsh .* nsh_grid{idim};
end
nsh = nsh(:)';

for idim = 1:msh.ndim
  csize = ones (1, 2*msh.ndim);
  csize([idim, msh.ndim+idim]) = [sp_univ(idim).nsh_max, msh.nel_dir(idim)];
  crep = [sp_univ.nsh_max, msh.nel_dir];
  crep([idim, msh.ndim+idim]) = 1;

  conn{idim} = reshape (sp_univ(idim).connectivity(:,elem_list{idim}), csize);
  conn{idim} = repmat (conn{idim}, crep);
  conn{idim} = reshape (conn{idim}, [], msh.nel);
end

connectivity = zeros (space.nsh_max, msh.nel);
indices = ones (size (conn{1}));
for idim = 1:msh.ndim
  indices = indices & conn{idim} ~= 0;
end
for idim = 1:msh.ndim
  conn{idim} = conn{idim}(indices);
end
connectivity(indices) = sub2ind ([space.ndof_dir, 1], conn{:}); % The extra 1 makes things work in any dimension
connectivity = reshape (connectivity, space.nsh_max, msh.nel);

clear conn csize crep indices

sp = struct('nsh_max', space.nsh_max, 'nsh', nsh, 'ndof', space.ndof,  ...
            'ndof_dir', space.ndof_dir, 'connectivity', connectivity, ...
            'ncomp', 1, 'degree', space.degree);

shp = cell(1,msh.ndim);
for idim = 1:msh.ndim
  ssize = ones (1, 3*msh.ndim);
  ssize([idim, msh.ndim+idim, 2*msh.ndim+idim]) = [msh.nqn_dir(idim), sp_univ(idim).nsh_max, msh.nel_dir(idim)];
  srep = [msh.nqn_dir, sp_univ.nsh_max, msh.nel_dir];
  srep([idim, msh.ndim+idim, 2*msh.ndim+idim]) = 1;
  shp{idim} = reshape (shape_functions_der(:,:,elem_list{idim}), ssize);
  shp{idim} = repmat (shp{idim}, srep);
  shp{idim} = reshape (shp{idim}, msh.nqn, space.nsh_max, msh.nel);  
end

  sp.shape_functions = 1;
  for idim = 1:msh.ndim
    sp.shape_functions = sp.shape_functions .* shp{idim};
  end

end