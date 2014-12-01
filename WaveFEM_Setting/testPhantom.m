function K =  testPhantom(node,ratio,type)

% [node, elem] = squaremesh([0,2,0,1],0.005);
% ratio = 2;
% type = 4;
% transform the coordinates to reference domain [0,ratio]*[0,1]
node(:,1) = (node(:,1) - min(node(:,1)))/(max(node(:,1))-min(node(:,1)))*ratio;
node(:,2) = (node(:,2) - min(node(:,2)))/(max(node(:,2))-min(node(:,2)));

nC = length(node);
K = 1500*ones(nC,1);

if type == 1 % line + 1/4 ellipse
   % first layer
   Ind = false(nC,1);
   Ind(node(:,2)>= (-0.2/ratio)*node(:,1) + 0.5) = true;
   K(Ind) = 1600;
    
   % second layer 
   Ind = false(nC,1);
   Ind((node(:,1)-ratio).^2*0.06 + (node(:,2)-1).^2 <= 0.25^2) = true;
   K(Ind) = 1900;
   
elseif type == 2  % two lines
    K = 2000*ones(nC,1);
   % first layer
   Ind = false(nC,1);
   Ind(node(:,2)<= 1/(4*ratio)*node(:,1) + 0.5) = true;
   K(Ind) = 1600;
    
   % second layer
   Ind = false(nC,1);
   Ind(node(:,2)<= -0.2/(ratio)*node(:,1) + 0.3) = true;
   K(Ind) = 1500;
   
elseif type == 3  % one line + one circle
   % first layer
    K = 1600*ones(nC,1);
   Ind = false(nC,1);
   Ind(node(:,2)<= 1/(4*ratio)*node(:,1) + 0.5) = true;
   K(Ind) = 1500;
   
   % Inclusion
   Ind = false(nC,1);
   Ind((node(:,1)-ratio*0.5).^2 + (node(:,2)-0.25).^2 <= 0.1^2) = true;
   K(Ind) = 1900;
   
elseif type == 4 % four circles as inclusions
   % Inclusion
   Ind = false(nC,1);
   Ind((node(:,1)-ratio*0.2).^2 + (node(:,2)-0.2).^2 <= 0.07^2) = true;
   K(Ind) = 1900;
   
   Ind = false(nC,1);
   Ind((node(:,1)-ratio*0.4).^2 + (node(:,2)-0.4).^2 <= 0.1^2) = true;
   K(Ind) = 1800;
   
   Ind = false(nC,1);
   Ind((node(:,1)-ratio*0.6).^2 + (node(:,2)-0.6).^2 <= 0.1^2) = true;
   K(Ind) = 1900;
   
    Ind = false(nC,1);
   Ind((node(:,1)-ratio*0.8).^2 + (node(:,2)-0.8).^2 <= 0.07^2) = true;
   K(Ind) = 1800;
   
elseif type == 5   % two ribs phantom
   % layer
   Ind = false(nC,1);
   Ind(node(:,2)<= 0.5 & node(:,2)>= 0.2) = true;
   K(Ind) = 1600;
    
   % two ribs
   
   Ind = false(nC,1);
   Ind((node(:,1)-ratio*0.3).^2/4 + (node(:,2)-0.8).^2*2 <= 0.09^2) = true;
   Ind((node(:,1)-ratio*0.7).^2/4 + (node(:,2)-0.8).^2*2 <= 0.09^2) = true;
   K(Ind) = 2500;
    
    
else
   K = 1500*ones(nC,1);
end

% showsolution(node,elem,K)
% view(0,-90)
end

