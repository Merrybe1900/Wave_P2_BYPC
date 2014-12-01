function test_phantom_small
xmin = -0.05;
xmax = 0.05;
ymin = -0.005;  % emitting b.c.
ymax = 0.035;

dx = 0.5000e-03;         % initial coarse meshsize [m]
[node,elem] = squaremesh([xmin,xmax,ymin,ymax], dx); % right triangle uniform mesh
soundSpeed = SoundSpeedPhantomBYsmall_ribs(node, 1);

showsolution(node,elem,soundSpeed);
view([0,-90])
caxis([1400,1600])
end