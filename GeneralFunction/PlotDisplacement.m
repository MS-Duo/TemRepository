function PlotDisplacement(Coordinate,Nodes,U)
figure;
hold on;

patch('Faces',Nodes(:,1:4),'Vertices',Coordinate,...
      'FaceVertexCData',U,...
      'FaceColor','interp','EdgeColor','none');

colormap(AbaqusColorMap());
axis equal;
xlabel('X'); ylabel('Y');

Co = colorbar('eastoutside');
set(gca, 'FontName','Times New Roman', 'FontSize',14)
Co.Label.String = 'U, Magnitude';

NTicks = 12;
[cmin, cmax] = caxis;
Co.Ticks = linspace(cmin,cmax,NTicks);
yt_labels = arrayfun(@(x) sprintf('%.3e', x), Co.Ticks, 'UniformOutput', false);
Co.YTickLabel = yt_labels;
end
