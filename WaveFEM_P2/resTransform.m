function [resT,E] = resTransform(uRef,u,p)

uRefH = transform(uRef);
uH = transform(u);
esyn = (uRefH.*uRefH + uRef.*uRef).^(1/2);
eobs = (uH.*uH + u.*u).^(1/2);

E = esyn.^p - eobs.^p;

EH = transform(E.*uRefH.*esyn.^(p-2));
resT = E.*uRef.*esyn.^(p-2) - EH;
end