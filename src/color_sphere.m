% Create a partial color-coded sphere for interpretation of fibre
% orientations.

background_color = 'k';

[X,Y,Z] = sphere(250);
idx = X<0 | Y<0 | Z<0;
X(idx)=NaN;Y(idx)=NaN;Z(idx)=NaN;
C=cat(3,X,Y,Z);

figure('Color',background_color)
hold on
surf(X,Y,Z,C,'EdgeColor','none')
quiver3( 1.02,0,0,0.2,0,0,'r','LineWidth',2,'ShowArrowHead','off')
quiver3( 0,1.02,0,0,0.2,0,'g','LineWidth',2,'ShowArrowHead','off')
quiver3( 0,0,1.02,0,0,0.2,'b','LineWidth',2,'ShowArrowHead','off')
text(1.3,0,0,'x','Color','r')
text(0,1.3,0,'y','Color','g')
text(0,0,1.3,'z','Color','b')
set(findobj(gca,'Type','Text'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',24,'FontWeight','bold',...
    'FontName','Arial')
axis equal tight off
set(gca,'Color',background_color)
view(120,25)