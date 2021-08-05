% Create a partial color-coded sphere for interpretation of fibre
% orientations.

background_color = 'w';
% xtext = 'x';
% ytext = 'y';
% ztext = 'z';
xtext = 'L-R';
ytext = 'A-P';
ztext = 'S-I';

sgnx = -1;
sgny = -1;

[X,Y,Z] = sphere(250);
if sgnx < 0;idx1 = X>0;else;idx1 = X<0;end
if sgny < 0;idx2 = Y>0;else;idx2 = Y<0;end
% idx = idx1 | idx2 | Z<0;
idx=[];
X(idx)=NaN;Y(idx)=NaN;Z(idx)=NaN;
C=cat(3,abs(X),abs(Y),abs(Z));

figure('Color',background_color)
hold on
surf(X,Y,Z,C,'EdgeColor','none')
quiver3( sgnx*1.02,0,0,sgnx*0.2,0,0,'r','LineWidth',2,'ShowArrowHead','off')
quiver3( 0,sgny*1.02,0,0,sgny*0.2,0,'g','LineWidth',2,'ShowArrowHead','off')
quiver3( 0,0,1.02,0,0,0.2,'b','LineWidth',2,'ShowArrowHead','off')
text(sgnx*1.3,0,0,xtext,'Color','r')
text(0,sgny*1.3,0,ytext,'Color','g')
text(0,0,1.3,ztext,'Color','b')
set(findobj(gca,'Type','Text'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',48,'FontWeight','bold',...
    'FontName','Arial')
axis equal tight off
set(gca,'Color',background_color)
% view(120,25)
view(-40, 10)
% view(140, 10)
material dull
set(gca,'clipping','off')