function fig=plot2dSolutionTimeSection(u,space,geometry,nrows,ncols,npts)

switch nargin
    case 4
        nrows=3;
        ncols=3;
        npts=50;
    case 5
        ncols=3;
        npts=50;
    case 6
        npts=50;
end

x={linspace(0,1,npts)};
t=linspace(0,0.75,nrows*ncols);
[eu, F] = sp_eval (u, space, geometry, [repmat(x,[1,numel(space.knots)-1]),{t}]);
min_eu=min(eu(:));
max_eu=max(eu(:));

fig=figure();
tiledlayout(nrows,ncols,'TileSpacing','compact','Padding','tight');
if numel(space.knots)==3
    for it=1:numel(t)
        nexttile;
        hold on;
        surf(squeeze(F(1,:,:,1)),squeeze(F(2,:,:,1)),eu(:,:,it));
        colormap(jet);
        fig.CurrentAxes.CLim=[min_eu,max_eu];
        fig.CurrentAxes.ZLim=[min_eu,max_eu];
        shading interp;
        view(2);
        axis tight;
        axis equal;
        title('',sprintf('t=%g',F(end,1,1,it)),'FontSize',8);
        box on;
    end
    cb = colorbar;
    cb.Layout.Tile = 'north';
elseif numel(space.knots)==2
    for it=1:numel(t)
        nexttile;
        hold on;
        plot(reshape(F(1,:,1),[],1),eu(:,it),'LineWidth',1.5);
        fig.CurrentAxes.YLim=[min_eu,max_eu];
        shading interp;
        title('',['t=',num2str(F(end,1,it))]);
        box on;
    end
end

end