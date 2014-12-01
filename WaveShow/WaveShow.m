function WaveShow
close all

dx = 2.0e-3;
xmin = -0.08;
xmax = 0.08;
ymin = -0.022;
ymax = 0.042;

maxRefine = 1;

recStep = 110;
source_sf = strcat('.\WaveRecData\reconstruction_P2_PE', num2str(recStep),'.mat');
save_sf_preflix = '.\WaveFig\reconstruction_PE2_';
save_fig = true;
save_eps = false;

x0 = -0.002;
y0 = 0.0035;

plotLimX = [1400,1640];
plotLimY = [1400,1640];

%*****************************************************************%
[node,elem] = squaremesh([xmin,xmax,ymin,ymax], dx); % right triangle uniform mesh
for i = 1 : maxRefine
  [node, elem] = uniformrefine(node,elem);
end
%-----------Split each p2 element into 6 small p1 elements ------%
bary = (node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:))/3;
[elem2dof,edge,~] = dofP2(elem); % convert to 6 points P2 elem
edgeMid = (node(edge(:,1),:) + node(edge(:,2),:))/2;
nodeP2= [node;edgeMid;bary]; % convert to masslumping nodes

nD = size(node,1) + size(edge,1);
nE = size(elem,1);
elem2dof = [elem2dof,(nD+1:nD+nE)']; % convert to mass lumping p2 elems

% split the p2 elements into p1 triangles
elemR = [elem2dof(:,1),elem2dof(:,6),elem2dof(:,7);...
         elem2dof(:,1),elem2dof(:,7),elem2dof(:,5);...
         elem2dof(:,5),elem2dof(:,7),elem2dof(:,3);...
         elem2dof(:,7),elem2dof(:,4),elem2dof(:,3);...
         elem2dof(:,7),elem2dof(:,2),elem2dof(:,4);...
         elem2dof(:,6),elem2dof(:,2),elem2dof(:,7)];
     
K = SoundSpeedPhantomBY_ribs(nodeP2,1);
% ---- show scan image ----- %
nC = size(node,1);
load(source_sf, 'f');
pic1 = figure(1);
showsolution(nodeP2,elemR, 1500./sqrt(1+f));
%showsolution(nodeP2,elemR,  f/0.121);
view([0,-90]);
caxis([1400,1600])
colormap gray
hold on

nodex = find(abs(nodeP2(:,1)- x0)<dx/(2^(maxRefine+5)));
if isempty(nodex)
    error('choose another x0');
end

[nodexs,nodex1] = sort(nodeP2(nodex,2));
nodex1 = nodex(nodex1);

plot(x0*ones(length(nodexs),1), nodexs,'r');
hold on;

nodey = find(abs(nodeP2(:,2) - y0)<dx/(2^(maxRefine+5)));
if isempty(nodey)
    error('choose another y0');
end
[nodeys,nodey1] = sort(nodeP2(nodey,1));
nodey1 = nodey(nodey1);

plot(nodeys,y0*ones(length(nodeys),1),'b');
hold off

if save_fig
    sf1 = strcat(save_sf_preflix, num2str(recStep),'.fig');
    saveas(pic1,sf1);
end
if save_eps
    sf1 = strcat(save_sf_preflix, num2str(recStep),'.eps');
    saveas(pic1,sf1,'epsc2');
end
% ------- Show vertical profile ----%
pic2 = figure(2);
axes('FontSize',16);
ylim(plotLimY);
plot(nodexs,1500./sqrt(1+f(nodex1)),'r','Linewidth',1.5);
hold on
plot(nodexs,K(nodex1),'k');
axis([-0.024,0.044,1420,1620]);
xlabel(strcat('y line,', 'x==',num2str(x0)));
ylabel('SoundSpeed');
grid on
hold off
if save_fig
   saveas(pic2, strcat(save_sf_preflix, num2str(recStep),'_ProfileX.fig'));
end
if save_eps
   saveas(pic2, strcat(save_sf_preflix, num2str(recStep),'_ProfileX.eps'),'epsc2');
end
% ------- Show horizontal profile ----%

pic3 =  figure(3);
axes('FontSize',16);
ylim(plotLimX);
plot(nodeys,1500./sqrt(1+f(nodey1)),'b','Linewidth',1.5);
hold on
plot(nodeys,K(nodey1),'k');
axis([-0.1,0.1,1420,1620]);
xlabel(strcat('x line,', 'y==',num2str(y0)));
ylabel('SoundSpeed');
grid on
hold off
if save_fig
   saveas(pic3, strcat(save_sf_preflix, num2str(recStep),'_ProfileY.fig'));
end
if save_eps
   saveas(pic3, strcat(save_sf_preflix, num2str(recStep),'_ProfileY.eps'),'epsc2');
end
% -------- Error Plot ----------%
% Offset plot

pic4 =  figure(4);
showsolution(nodeP2,elemR, abs(1500./sqrt(1+f)-K)./K);
view([0,-90]);
caxis([0,0.04]);
colorbar

if save_fig
   saveas(pic4, strcat(save_sf_preflix, num2str(recStep),'_ErrorMap.fig'));
end

if save_eps
   saveas(pic4, strcat(save_sf_preflix, num2str(recStep),'_ErrorMap.eps'),'epsc2');
end
end