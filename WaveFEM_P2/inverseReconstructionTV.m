function [f,h]= inverseReconstructionTV(mesh,bda,source,signal,A,As,area,f,uSOR,Nt,...
             dt,displayOption,count,h,sf_Mesure2)
%**************** reconstruction function ***************%
% Author: Yun Bai
% Date: 21.08.2014
%*********************************************************%

% assemble mass matrices for global and boundary
Kn = mesh.K0^2*(1+f);
Ke = 1/3*(Kn(mesh.elem(:,1)) + Kn(mesh.elem(:,2)) + Kn(mesh.elem(:,3)));
[M,Me,Ms] = massMatrixWaveP2(mesh,bda,area,Ke,Kn,dt);
beta = 0.2;

uhS = cell(signal.numS,1);
residual = zeros(mesh.Ndof,1);
mesh.solType = 'i';

if mesh.quadorder == 3.5
    
       parfor ks = 1 : signal.numS
              % forward solver %
              
              [uhSOR,uDir,u0,u1] = waveForwardSolverP2massLumping2(mesh,bda,source,signal,...
                                    ks,A,As,M,Me,Ms,Nt,dt,displayOption,Kn);
              uhS{ks} = uSOR{ks} - uhSOR;
    %         residual = residual  +  waveBackwardSolverP2massLumping(mesh,bda, uhS{ks},uDir,...
    %                       u,u0,A,As,M,Me,Ms,Nt,dt,displayOption); 
              residual = residual  +  waveBackwardSolverP2massLumping2(mesh,bda, uhS{ks},uDir,...
                                    u0,u1,A,As,M,Me,Ms,Nt,dt,displayOption); 
       end
       TV = zeros(mesh.Ndof,1);
       TV(bda.freenode) =...
       A(bda.freenode,bda.freenode)*f(bda.freenode)/sqrt(f(bda.freenode)'*A(bda.freenode,bda.freenode)*f(bda.freenode)+ beta);
       
       residual = residual + 500*TV;
       
else %
       parfor ks = 1 : signal.numS   
           % forward solver %
           [uhSOR,uDir,u,u0] = waveForwardSolverP2(mesh,bda,source,signal,...
                                    ks,A,As,M,Me,Ms,Nt,dt,displayOption,Kn);
           % Backward solver %
           uhS{ks} = uSOR{ks} - uhSOR;
           residual = residual  +  waveBackwardSolverP2(mesh,bda, uhS{ks},uDir,...
                                     u,u0,A,As,M,Me,Ms,Nt,dt,displayOption); 
       end

 end

%%  line search %%
if ~mod(count,5)|| count == 1
    disp('-- line search --');

    alpha1 = 500;
    alpha3 = 0.3;
    
    mR = max(abs(residual));
    if  mR >=  1e-3
        h = alpha1*(mesh.K0)^2;
    else
        showsolution(mesh.node, mesh.elem, 1500./sqrt(1+f(1:mesh.nC)));
        view([0,-90]);
        pause(0.05);
        error('residual is too small');
    end
    
    %r0 = sum(residual);

    uhSe = zeros(bda.lenSOR-1,Nt);
    for i =1:signal.numS
        uhSe = uhSe + 1/6*(uhS{i}(1:bda.lenSOR-1,:).^2 ...
                    + 4*uhS{i}(bda.lenSOR+1:2*bda.lenSOR-1,:).^2 ...
                    + uhS{i}(2:bda.lenSOR,:).^2);
    end

    costF0 = sum(bda.elSOR'*uhSe)*(dt*(1/mesh.K0)^2); 
    c1 =costF0*0.000001;

    mesh.maxIterL = 3;

    for iter = 1:mesh.maxIterL
         disp(strcat('line search iteration: ',num2str(iter)));
         ft = f;
         ft  = ft + h.*residual;
         Ktn = mesh.K0^2*(1 + ft);
         Kte = 1/3*(Ktn(mesh.elem(:,1)) + Ktn(mesh.elem(:,2)) + Ktn(mesh.elem(:,3)));

         [M,Me,Ms] = massMatrixWaveP2(mesh,bda,area,Kte,Ktn,dt);
         mesh.solType = 'e';
         if mesh.quadorder == 3.5
              As = stiffMatrixABCP2(mesh,bda,Ktn);
             parfor ks =1:signal.numS
                    uhSOR = waveForwardSolverP2massLumping2(mesh,bda,source,signal,...
                                        ks,A,As,M,Me,Ms,Nt,dt,displayOption,Ktn);
                    uhS{ks} = uSOR{ks} - uhSOR;

             end
         else
             parfor ks =1:signal.numS
                    uhSOR = waveForwardSolverP2(mesh,bda,source,signal,...
                                        ks,A,As,M,Me,Ms,Nt,dt,displayOption,Ktn);
                    uhS{ks} = uSOR{ks} - uhSOR;

             end
         end

        uhSe = zeros(bda.lenSOR-1,Nt);

        for i =1:signal.numS
            uhSe = uhSe + 1/6*(uhS{i}(1:bda.lenSOR-1,:).^2 ...
                        + 4*uhS{i}(bda.lenSOR+1:2*bda.lenSOR-1,:).^2 ... 
                        + uhS{i}(2:bda.lenSOR,:).^2);
        end
        costFh = sum(bda.elSOR'*uhSe)*(dt*(1/mesh.K0)^2); 
        err = costFh - costF0 + c1*sum(h.*residual);
        if err <= 0
           break;
        else
           h = alpha3*h;
        end
    end
end
    
f = f + h.*residual;

% impose a hard upper bound 
Kn = 1./(mesh.K0*sqrt((1+f)));
indicator = false(mesh.Ndof,1);
indicator(Kn >= 2700) = true;
f(indicator) = -0.6672;
indicator = false(mesh.Ndof,1);
indicator(Kn <= 1450) = true;
f(indicator) =  0.0702;

end

