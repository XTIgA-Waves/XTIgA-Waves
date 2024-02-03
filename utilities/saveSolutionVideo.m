function fig=saveSolutionVideo(u,space,geometry,fileName,nframe,npts)

switch nargin
    case 4
        nframe=10;
        npts=50;
    case 5
        npts=50;
end

x={linspace(0,1,npts)};
t=linspace(0,1,nframe);
if geometry.rdim==3
    [eu, F] = sp_eval (u, space, geometry, [x,x,{t}]);
else
    [eu, F] = sp_eval (u, space, geometry, [x,{t}]);
end
min_eu=min(eu(:));
max_eu=max(eu(:));

fileName=sprintf('%s.avi',fileName);
video = VideoWriter(fileName);
video.Quality = 100;
video.FrameRate = 10;
open(video);
fig=figure('units','pixels','position',[200 100 720 720]);

colormap(jet);
for iframe=1:numel(t)
    tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
    if geometry.rdim==3
        nexttile;
        surf(squeeze(F(1,:,:,1)),squeeze(F(2,:,:,1)),eu(:,:,iframe));
        fig.CurrentAxes.CLim=[min_eu,max_eu];
        fig.CurrentAxes.ZLim=[min_eu,max_eu];
        shading interp;
        light;
        nexttile;
        surf(squeeze(F(1,:,:,1)),squeeze(F(2,:,:,1)),eu(:,:,iframe));
        view(2);
        shading interp;
        sgtitle(sprintf('t = %.4g',F(end,1,1,it)),'FontSize',28);
        fig.CurrentAxes.CLim=[min_eu,max_eu];
        fig.CurrentAxes.ZLim=[min_eu,max_eu];
        cb = colorbar;
        cb.Layout.Tile = 'south';
    else
        plot(squeeze(F(1,:,:,1)),eu(:,iframe));
        sgtitle(sprintf('t = %.4g',F(end,1,iframe)),'FontSize',28);
        fig.CurrentAxes.YLim=[min_eu,max_eu];
    end

    axis square;


    writeVideo(video,getframe(gcf));
end
close(video);

end