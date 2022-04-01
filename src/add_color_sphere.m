function add_color_sphere(O,scale)
% Create a partial color-coded sphere for interpretation of fibre
% orientations.

% background_color = 'k';
fontsize = 16;
[X,Y,Z] = sphere(250);
% idx = X<0 | Y<0 | Z<0;
% X(idx)=NaN;Y(idx)=NaN;Z(idx)=NaN;
C=cat(3,X,Y,Z);
X=X*scale+O(1);
Y=Y*scale+O(2);
Z=Z*scale+O(3);

% figure('Color',background_color)
hold on
surf(X,Y,Z,C,'EdgeColor','none')
quiver3( O(1)+scale,O(2),O(3),0.2*scale,0,0,'r','LineWidth',2,'ShowArrowHead','off')
quiver3( O(1),O(2)+scale,O(3),0,0.2*scale,0,'g','LineWidth',2,'ShowArrowHead','off')
quiver3( O(1),O(2),O(3)+scale,0,0,0.2*scale,'b','LineWidth',2,'ShowArrowHead','off')
text(O(1)+1.3*scale,O(2),O(3),'x','Color','r')
text(O(1),O(2)+1.3*scale,O(3),'y','Color','g')
text(O(1),O(2),O(3)+1.3*scale,'z','Color','b')
set(findobj(gca,'Type','Text'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',fontsize,'FontWeight','bold',...
    'FontName','Arial')
axis equal tight off
% set(gca,'Color',background_color)
% view(120,25)
material dull