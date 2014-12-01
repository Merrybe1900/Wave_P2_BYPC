function [M,b] = P2Basis_test
%--------- boundary basis function test -----------------%
x = [0,1,0.5]';

lambdagR = [1-x,x];
bdphi(:,1) = (2*lambdagR(:,1)-1).*lambdagR(:,1);
bdphi(:,2) = (2*lambdagR(:,2)-1).*lambdagR(:,2);
bdphi(:,3) = 4*lambdagR(:,1).*lambdagR(:,2);

Md = zeros(3,3);
weight = [1/6,1/6,2/3]';
for i = 1:3
    for j = 1:3
        pij = sum(weight.*bdphi(:,i).*bdphi(:,j));
        Md(i,j) = pij;
    end
end

%--------- basis function test -----------------%
x = [1,0,0,0,0.5,0.5,1/3]';
y = [0,1,0,0.5,0,0.5,1/3]';

lambda3 = 1-x-y;
lambda1 = x;
lambda2 = y;

dlam3 = [-ones(7,1),-ones(7,1)];
dlam1 = [ones(7,1),zeros(7,1)];
dlam2 = [zeros(7,1),ones(7,1)];

weight = [1/20,1/20,1/20,2/15,2/15,2/15,9/20]';

p = zeros(length(weight),7);
p(:,1) = 2*lambda1.^2-lambda1 + 3.*lambda1.*lambda2.*lambda3;
p(:,2) = 2*lambda2.^2-lambda2 + 3.*lambda1.*lambda2.*lambda3;
p(:,3) = 2*lambda3.^2-lambda3 + 3.*lambda1.*lambda2.*lambda3;

p(:,4) = 4.*lambda2.*lambda3 - 12.*lambda1.*lambda2.*lambda3;
p(:,5) = 4*lambda1.*lambda3 - 12*lambda1.*lambda2.*lambda3;
p(:,6) = 4*lambda1.*lambda2 - 12*lambda1.*lambda2.*lambda3;
p(:,7) = 27*lambda1.*lambda2.*lambda3;

M = zeros(7,7);
for i = 1:7
    for j = 1:7
        pij = sum(weight.*p(:,i).*p(:,j))/2;
        M(i,j) = pij;
    end
end

b= zeros(7,1);
for i = 1:7
    b(i) = sum(weight.*p(:,i))/2;
end

dp = zeros(length(weight),2,7);
lambda1 = repmat(lambda1,1,2);
lambda2 = repmat(lambda2,1,2);
lambda3 = repmat(lambda3,1,2);
dp(:,:,1) = 4*dlam1.*lambda1- dlam1 + 3.*(dlam1.*lambda2.*lambda3 ...
          + dlam2.*lambda1.*lambda3 + dlam3.*lambda1.*lambda2);
dp(:,:,2) = 4*dlam2.*lambda2- dlam2 + 3.*(dlam1.*lambda2.*lambda3 ...
          + dlam2.*lambda1.*lambda3 + dlam3.*lambda1.*lambda2);
dp(:,:,3) = 4*dlam3.*lambda3- dlam3 + 3.*(dlam1.*lambda2.*lambda3 ...
          + dlam2.*lambda1.*lambda3 + dlam3.*lambda1.*lambda2);


dp(:,:,4) = 4*(dlam2.*lambda3 + dlam3.*lambda2) -  12.*(dlam1.*lambda2.*lambda3 ...
          + dlam2.*lambda1.*lambda3 + dlam3.*lambda1.*lambda2);
dp(:,:,5) = 4*(dlam1.*lambda3 + dlam3.*lambda1) -  12.*(dlam1.*lambda2.*lambda3 ...
          + dlam2.*lambda1.*lambda3 + dlam3.*lambda1.*lambda2);
dp(:,:,6) = 4*(dlam1.*lambda2 + dlam2.*lambda1) -  12.*(dlam1.*lambda2.*lambda3 ...
          + dlam2.*lambda1.*lambda3 + dlam3.*lambda1.*lambda2);
dp(:,:,7) = 27.*(dlam1.*lambda2.*lambda3 ...
          + dlam2.*lambda1.*lambda3 + dlam3.*lambda1.*lambda2);
      
A = zeros(7,7);
weight = repmat(weight,1,2);
for i = 1:7
    for j = 1:7
        pij = sum(sum(dp(:,:,j).*(weight.*dp(:,:,i))))/2;
        A(i,j) = pij;
    end
end
end