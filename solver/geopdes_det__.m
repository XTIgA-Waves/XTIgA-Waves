% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2015 Rafael Vazquez
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

function d = geopdes_det__ (v)

vsize = size (v);

if (vsize(1) == 1 && vsize(2) == 1)
    d = v(1,1,:,:);

elseif (vsize(1) == 2 && vsize(2) == 2)
    d = v(1,1,:,:) .* v(2,2,:,:) - v(2,1,:,:) .* v(1,2,:,:);

elseif (vsize(1) == 3 && vsize(2) == 3)
    d = v(1,1,:,:) .* (v(2,2,:,:) .* v(3,3,:,:) - v(2,3,:,:) .* v(3,2,:,:))...
        + v(1,2,:,:) .* (v(2,3,:,:) .* v(3,1,:,:) - v(2,1,:,:) .* v(3,3,:,:))...
        + v(1,3,:,:) .* (v(2,1,:,:) .* v(3,2,:,:) - v(2,2,:,:) .* v(3,1,:,:));

elseif (vsize(1) == 4 && vsize(2) == 4)
    sigma=-1;
    Sn=perms(1:4);
    d=zeros([1,1,prod(size(v,3:ndims(v)))]);
    for iperm=size(Sn,1):-1:1
        sigma=-sigma;
        d=d+sigma*v(1,Sn(iperm,1),:) .*v(2,Sn(iperm,2),:) .*v(3,Sn(iperm,3),:) .*v(4,Sn(iperm,4),:);
    end
    d=reshape(d,size(v,3:ndims(v)));
    
elseif (vsize(1) > vsize(2))
    % G = v^t * v, first fundamental form for v = DF
    for ii = 1:vsize(2)
        for jj = 1:vsize(2)
            G(ii,jj,:,:) = sum (v(:,ii,:,:).*v(:,jj,:,:), 1);
        end
    end

    d = sqrt (geopdes_det__ (G));
end

d = squeeze (d);

end
