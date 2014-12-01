function [resT,E] = resTransformLog(uRef,u)

uRefH = transform(uRef);
uH = transform(u);
esyn = (uRefH.*uRefH + uRef.*uRef).^(1/2);
eobs = (uH.*uH + u.*u).^(1/2);

E = log(esyn) - log(eobs);

EH = transform(E.*uRefH.*esyn.^(-2));
resT = - E.*uRef.*esyn.^(-2)+ EH;
end