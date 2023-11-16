function [errh1, errl2, errh1s, errh1_elem, errl2_elem, errh1s_elem] = sp_h1_c_error_tp (space, msh, u, uex, graduex, c)

if (numel(u) ~= space.ndof)
    error ('Wrong size of the vector of degrees of freedom')
end

errl2 = 0; errh1s = 0;
errh1_elem = zeros (1, msh.nel); errl2_elem = zeros (1, msh.nel); errh1s_elem = zeros (1, msh.nel);

for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', true, 'gradient', true);

    [~, err_l2, err_h1s, errh1_col_elem, errl2_col_elem, errh1s_col_elem] = sp_h1_c_error (sp_col, msh_col, u, uex, graduex, c);

    errh1s = errh1s + err_h1s.^2;
    errl2 = errl2 + err_l2.^2;

    errh1_elem(:,msh_col.elem_list)  = errh1_col_elem;
    errl2_elem(:,msh_col.elem_list)  = errl2_col_elem;
    errh1s_elem(:,msh_col.elem_list) = errh1s_col_elem;
end

errh1 = sqrt (errl2 + errh1s);
errl2 = sqrt (errl2);
errh1s = sqrt (errh1s);

end
