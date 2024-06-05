function [errh1, errl2, errh1s] = sp_h1_error_spline_ex (space, msh, u, space_fine, geometry, u_fine,varargin)

dim=geometry.rdim;
errh1sDirSquared=0;
errl2sDirSquared=0;

parfor iDir1=1:numel(msh.qn{1})
    if dim==2
        jacobian=feval(msh.map_der,{msh.qn{1}(iDir1),msh.qn{2}(:)});
        [grad_valex,F] = sp_eval(u_fine,space_fine,geometry,{msh.qn{1}(iDir1),msh.qn{2}},{'value','gradient'});
        grad_valu = sp_eval(u,space,geometry,{msh.qn{1}(iDir1),msh.qn{2}},{'value','gradient'});
    elseif dim==3
        jacobian=msh.map_der({msh.qn{1}(iDir1),msh.qn{2},msh.qn{3}});
        [grad_valex,F] = sp_eval(u_fine,space_fine,geometry,{msh.qn{1}(iDir1),msh.qn{2},msh.qn{3}},{'value','gradient'});
        grad_valu = sp_eval(u,space,geometry,{msh.qn{1}(iDir1),msh.qn{2},msh.qn{3}},{'value','gradient'});
    end
    jacDet=zeros(size(jacobian,3),1);
    for im=1:size(jacobian,3)
        jacDet(im)=det(jacobian(:,:,im));
    end

    qw = msh.qw{1}(iDir1);
    for idim = 2:msh.ndim
        qw = kron (msh.qw{idim}, qw);
    end

    jacDet=reshape(jacDet,size(qw));

    w = abs(qw.*jacDet);

    grad_valex{1} = reshape (grad_valex{1}, 1, 1, prod(msh.nqn_dir(2:end)), prod(msh.nel_dir(2:end)));
    grad_valex{2} = reshape (grad_valex{2}, space.ncomp, msh.rdim, prod(msh.nqn_dir(2:end)), prod(msh.nel_dir(2:end)));
    grad_valu{1} = reshape (grad_valu{1}, 1, 1, prod(msh.nqn_dir(2:end)), prod(msh.nel_dir(2:end)));
    grad_valu{2} = reshape (grad_valu{2}, space.ncomp, msh.rdim, prod(msh.nqn_dir(2:end)), prod(msh.nel_dir(2:end)));

    if numel(varargin)==1
        F=reshape(F,geometry.rdim,[]);
        x=[];
        for idim=1:geometry.rdim
            x{idim}=F(idim,:);
        end
        cweight=reshape(feval(varargin{1},x{:}),1,[]);
        cweight=repmat(cweight,[geometry.rdim-1,1]);
        cweight=cat(1,cweight,ones(1,size(cweight,2)));
        tmp1=reshape(vecnorm(reshape(grad_valu{2} - grad_valex{2},geometry.rdim,[]).*cweight,2,1),size(w)).^2.*w;
    else
        tmp1=reshape(vecnorm(reshape(grad_valu{2} - grad_valex{2},geometry.rdim,[]),2,1),size(w)).^2.*w;
    end
    errh1sDirSquared(iDir1)=sum(tmp1,'all');

    tmp2=squeeze(grad_valu{1} - grad_valex{1}).^2.*w;
    errl2sDirSquared(iDir1)=sum(tmp2,'all');
end

errh1sDirSquared=sum(errh1sDirSquared);
errl2sDirSquared=sum(errl2sDirSquared);

errh1s = sqrt (errh1sDirSquared);
errl2 = sqrt (errl2sDirSquared);
errh1 = errl2 + errh1s;

end