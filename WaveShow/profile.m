function profile(node,x0,y0,f)
nodex = find(node(:,1) == x0);
nodey = find(node(:,1) == y0);

f0 = SoundSpeedPhantomBY_ribs(node,1);

figure(1)
plot(node(nodex,2),1500./sqrt(f(nodex)),'b');
hold on
plot(node(nodex,2),1500./sqrt(f0(nodex)),'k');
hold off

figure(2)
plot(node(nodey,1),1500./sqrt(f(nodey)),'r');
hold on
plot(node(nodey,1),1500./sqrt(f0(nodey)),'k');
hold off


end