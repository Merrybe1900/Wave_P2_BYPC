function soundSpeed = SoundSpeedPhantomBYsmall_ribs(node, withRibs)
if nargin < 2, withRibs = 0; end

%Sound speeds.
soundSpeedFat = 1470;                                                                                   % Fat sound speeds in [m/s]
soundSpeedTumor = 1600;                                                                                 % Tumor sound speeds in [m/s]
glandularTissueSoundSpeed = 1545;                                                                       % Grandular tissue sound speeds in [m/s]
soundSpeedMuscle = 1545;                                                                                % Muscle sound speeds in [m/s]
soundSpeedWater = 1500;                                                                                 % Water sound speeds in [m/s]
roi = struct('xMin', -0.02, 'xMax', 0.02, 'yMin', -0.05, 'yMax', 0.05, 'zMin', -0.005, 'zMax', 0.035);  % ROI is [xMin, xMax, yMin, yMax, zMin, zMax], e.g.  [-0.04, 0.04, -0.08, 0.08, -0.022, 0.042].

averageFilter = 1/4*ones(2,2);                                                                          % To average the phantom after generation.
homogenityLevel = 28;                                                                                   % In dB: Add Inhomogenity to the phantom by adding white Gaussian noise to the phantom.
radius = 4.5e-2;                                                                                          % Breast shape is a cut cylinder. Radius of the cylinder in [m].
centerY = 0;                                                                                            % Breast shape is a cut cylinder. Center of the cylinder in [m]
centerZ = 4.5e-2;                                                                                         % Breast shape is a cut cylinder. Center of the cylinder in [m]
roiCenter = [(roi.xMax+roi.xMin)/2 (roi.yMax+roi.yMin)/2 (roi.zMax+roi.zMin)/2];

%% Inclusion parameters
% Inclusion 1: Fat sphere position and sound speed.
x1 = 0;
y1 = -0.022;
z1 = 1.6e-2;
r1x = 3e-3;
r1y = 1.4e-3;
r1z = 1.0e-3;
soundSpeed1 = soundSpeedFat;

% Inclusion 2: Fat sphere position and sound speed.
x2 = -2e-2;
y2 = 1e-2;
z2 = -3e-2;
r2x = 0.35e-2;
r2y = 0.45e-2;
r2z = 0.35e-2;
soundSpeed2 = soundSpeedFat;

% Inclusion 3: Tumor 2mm close to 1mm tumors.
x3 = 0.0005;
y3 = -0.9e-2;
z3 = 1.7e-2;
r3x = 1.45e-3;
r3y = 1.45e-3;
r3z = 1.45e-3;
soundSpeed3 = soundSpeedTumor;

% Inclusion 4: Fat inclusion big one on bottom.
x4 = 0.5e-3;
y4 = 1.5e-2;
z4 = 1.5e-2;
r4x = 2e-3;
r4y = 3.5e-3;
r4z = 3.5e-3;
soundSpeed4 = soundSpeedFat;

% Inclusion 5: Tumor.
x5 = 0.5e-3;
y5 = 2.55e-2;
z5 = 1.605e-2;
r5x = 2.0e-3;
r5y = 2.0e-3;
r5z = 2.0e-3;
soundSpeed5 = soundSpeedTumor;


% Inclusion 7: Tumor.
x7 = 0;
y7 = -1.5e-2;
z7 = 1.2e-2;
r7x = 2e-3;
r7y = 1.5e-3;
r7z = 2.25e-3;
soundSpeed7 = soundSpeedTumor;

% Inclusion 8: Tumor 1mm size.
x8 = 0.000;
y8 = 0.006;
z8 = 0.014;
r8x = 0.75e-3;
r8y = 0.75e-3;
r8z = 0.75e-3;
soundSpeed8 = soundSpeedTumor;

% Inclusion 9: Tumor.
x9 = 0.9e-2;
y9 = 2.5e-2;
z9 = 2.1e-2;
r9x = 2.2e-3;
r9y = 3e-3;
r9z = 2.5e-3;
soundSpeed9 = soundSpeedTumor;

% Inclusion 10: Fat inclusion.
x10 = -1.0e-2;
y10 = -0.8e-2;
z10 = 1.3e-2;
r10x = 1.9e-3;
r10y = 2.1e-3;
r10z = 2.8e-3;
soundSpeed10 = soundSpeedFat;

% Inclusion 11: Tumor 1mm sphere.
x11 = 0.000;
y11 = 0.006;
z11 = 0.008;
r11x = 0.75e-3;
r11y = 0.75e-3;
r11z = 0.75e-3;
soundSpeed11 = soundSpeedTumor;

% Ribs
ribProfile = 1;                                                             % rib profile: 0 - rectangle; 1 - ellipse
soundSpeedRibs = 2500;                                                                % set the sound speed in ribs
ribDimensions = [(roi.xMax-roi.xMin)- 0.012 0.009 0.0038];                    % define ribs dimensions
nRibs = ceil(floor(abs(roi.yMax-roi.yMin-0.005) / ribDimensions(2)) / 2);   % let the ribs span as much as possible the y side of the ROI
ribExtensionAlongZ = 0.03 + [0 ribDimensions(3)];                                        % set the ribs' extent along the z coordinate
ribExtensionAlongX = ((roi.xMax-roi.xMin) - ribDimensions(1))/2 .* [1 -1] + ...
  [roi.xMin roi.xMax];                                                      % set the ribs' extent along the x coordinate
ribStartY = roiCenter(2) - (2*nRibs-1)/2*ribDimensions(2);                  % set the ribs' y extent (multiple ribs now)
for i = 1:nRibs
  ribExtensionAlongY(i,:) = ribStartY + 2 * (i-1) * ribDimensions(2) + ...
    [0 ribDimensions(2)];
end
numberOfPointsInsideRibs = 0;                                                              % count the number of points inside of ribs


%% Generate irregular boundary between the fat and grandular tissue.
% Each row has the format [centerx, centery, radius]
n = double(28);
theta = [0+(0:n-2)*(2*pi-0)/(floor(n)-1) 2*pi]+0.2312412;
bndCirc = zeros(length(theta), 3);
for i = 1 : length(theta)
    bndCirc(i,1) = centerY + 4.5e-2*cos(theta(i));
    bndCirc(i,2) = centerZ + 4.5e-2*sin(theta(i));
    bndCirc(i,3) = 0.45e-2;
end

%% Generate boundaries of the irregular tumor shape.
% Each row has the format [centerx, centery, radius]
n = double(6);
rStar = 4.625e-3;
theta = [0+(0:n-2)*(2*pi-0)/(floor(n)-1) 2*pi]+0.7312412;
cancerCirc = zeros(length(theta), 3);
for i = 1 : length(theta)
    cancerCirc(i,1) = 1.001e-2 + rStar *sin(theta(i));
    cancerCirc(i,2) = -0.3005e-2 + rStar *cos(theta(i));
    cancerCirc(i,3) = 2e-3;
end

%% Generate irregular shape boundary between the fat and the grandular tissue
n = double(6);
fatCirc = zeros(n, 3);
fatY = [-0.029 0.044];
for i = 1 : n
    fatCirc(i,1) = fatY(1) + (fatY(2)-fatY(1))/n*(i-1); % y
    fatCirc(i,2) = 0.02; % z
    fatCirc(i,3) = 5e-3; % radius
end

%% Muscle layers on the bottom.
muscle_a1 =     0.01963;
muscle_b1 =     0.105;
muscle_c1 =     0.06162;
muscle_a2 =           0;
muscle_b2 =    -0.01908;
muscle_c2 =   0.0003059;
muscle_a3 =     0.03013;
muscle_b3 =     -0.0824;
muscle_c3 =      0.2052;

%% For each point in the ROI check if it belongs to the previously defined structures.
nC = size(node,1);
soundSpeed = zeros(nC,1);
for p = 1:nC

    if size(node,2) == 2
        x = 0; 
        y = node(p,1);
        z = node(p,2);
    else
        x = node(p,1);
        y = node(p,2);
        z = node(p,3);
    end

    %% Set the flags to 1
    inside_flag = 1;
    fat_flag = 1;

    %% If it is below the breast phantom
    if z > 0.02
        % Check if it is in the muscle layer
        if (muscle_a1*exp(-((y-muscle_b1)/muscle_c1)^2) + muscle_a2*exp(-((y-muscle_b2)/muscle_c2)^2) + muscle_a3*exp(-((y-muscle_b3)/muscle_c3)^2) - z < 0)
            soundSpeed(p) = soundSpeedMuscle;
            inside_flag = 0;
        else  % Check if it in the fat layer
            for j = 1 : size(fatCirc,1)
                if (y-fatCirc(j,1))^2+(z-fatCirc(j,2))^2 <= fatCirc(j,3)^2
                    fat_flag  = 0;
                    break;
                end
            end
            if fat_flag
                soundSpeed(p) = soundSpeedFat;
                inside_flag = 0;
            else
                inside_flag = 1;
            end
        end
    else
        fat_flag = 0;
    end

    %% Check if it is outside the circles that defines the boundary between the fat and the grandular tissue.
    for j = 1 : size(bndCirc,1)
        if (y-bndCirc(j,1))^2+(z-bndCirc(j,2))^2 <= bndCirc(j,3)^2
            inside_flag = 0;
            break;
        end
    end

    %% Then it is inside the grandular tissue.
    if inside_flag == 1 && (y-centerY)^2 + (z-centerZ)^2 <= (radius)^2
        soundSpeed(p) = glandularTissueSoundSpeed;                                  % Glandular tissue
        if (z-z1)^2/r1z^2+(y-y1)^2/r1y^2 <= 2 && (x-x1)^2 < r1x^2
            soundSpeed(p) = soundSpeed1;
            %                 if (z-z2)^2/r2z^2+(y-y2)^2/r2y^2 <= 1 && (x-x2)^2 < r2x^2
            %                     soundSpeed(p) = soundSpeed2; %fat sphere
        elseif (z-z3)^2/r3z^2+(y-y3)^2/r3y^2 <= 1 && (x-x3)^2 < r3x^2
            soundSpeed(p) = soundSpeed3;                                            % Inclusion 3.
        elseif (z-z4)^2/r4z^2+(y-y4)^2/r4y^2 <= 1 &&  + (x-x4)^2 < r4x^2
            soundSpeed(p) = soundSpeed4;                                            % Inclusion 4.
        elseif (z-z5)^2/r5z^2+(y-y5)^2/r5y^2 <= 1 && (x-x5)^2 < r5x^2
            soundSpeed(p) = soundSpeed5;                                            % Inclusion 5.
        elseif (z-z7)^2/r7z^2+(y-y7)^2/r7y^2 <= 1 && (x-x7)^2 < r7x^2
            soundSpeed(p) = soundSpeed7;                                            % Inclusion 7.
        elseif (z-z8)^2/r8z^2+(y-y8)^2/r8y^2 < 1 && (x-x8)^2 < r8x^2
            soundSpeed(p) = soundSpeed8;                                            % Inclusion 8
        elseif (z-z9)^2/r9x^2+(y-y9)^2/r9y^2 <= 1 && (x-x9)^2 < r9x^2
            soundSpeed(p) = soundSpeed9;                                            % Inclusion 9.
        elseif (z-z10)^2/r10x^2+(y-y10)^2/r10y^2 <= 1 && (x-x10)^2 < r10x^2
            soundSpeed(p) = soundSpeed10;                                           % Inclusion 10.
        elseif (z-z11)^2/r11z^2+(y-y11)^2/r11y^2 + (x-x11)^2/r11x^2 < 1
            soundSpeed(p) = soundSpeed11;                                           % Inclusion 11.
        else                                                                        % Check the iregular tumor
            inside_tumor = 1;
            for j = 1 : size(cancerCirc,1)
                if (z-cancerCirc(j,1))^2+(y-cancerCirc(j,2))^2 <= cancerCirc(j,3)^2
                    inside_tumor = 0;
                    break;
                end
            end
            if inside_tumor == 1 && (z-1e-2)^2+(y+0.3e-2)^2 <= (rStar )^2 && (x - 0)^2 <= (3e-3)^2
                soundSpeed(p) = soundSpeed5;                                        % Inclusion 5.
            end
        end
    elseif ~fat_flag && (y-centerY)^2 + (z-centerZ)^2 > (radius+1.5e-3)^2           % Check if the point is outside all
        soundSpeed(p) = soundSpeedWater;
    elseif ~fat_flag                                                                % Now the point shoud be between the regular and irregular boundaries
        soundSpeed(p) = soundSpeed2;                                                % Subcutaneous fat
    end

    % set the sound speed to rib sound speed if inside any of the ribs
    insideRib = 0;

    if ribProfile == 1 % elliptic profile
        % check each ellipsoidal yz profile
        % any is good
        for n = 1:nRibs
            % determine ellipsoid parameters first
            yCenter = 1/2 * sum(ribExtensionAlongY(n,:));
            zCenter = 1/2 * sum(ribExtensionAlongZ);
            rY = abs(ribExtensionAlongY(n,2) - yCenter);
            rZ = abs(ribExtensionAlongZ(2) - zCenter);
            % check if in the ellipsoid
            if (y-yCenter).^2 ./ rY.^2 + (z-zCenter).^2 / rZ.^2 <= 1
                insideRib = 1;
                break;
            end
        end
        % check x limits
        if x < ribExtensionAlongX(1) || x > ribExtensionAlongX(2)
            insideRib = 0;
        end
    else % rectangular profile
        % check y limits
        for n = 1:nRibs
            % any is good
            if y >= ribExtensionAlongY(n,1) && y <= ribExtensionAlongY(n,2)
                insideRib = 1;
                break;
            end
        end
        % check x limits
        if x < ribExtensionAlongX(1) || x > ribExtensionAlongX(2)
            insideRib = 0;
        end
        % check z limits
        if z < ribExtensionAlongZ(1) || z > ribExtensionAlongZ(2)
            insideRib = 0;
        end
    end
    % set the sound speed if all
    if insideRib
        soundSpeed(p) = soundSpeedRibs;
        numberOfPointsInsideRibs = numberOfPointsInsideRibs + 1;
    end

end


%% Save the indexes of the water cells.
waterIndex = soundSpeed == soundSpeedWater;
ribsIndex = soundSpeed == soundSpeedRibs;
noRibsIndex = soundSpeed ~= soundSpeedRibs;

%% Add noise to all points.
% adding noise to the phantom
RandStream.setGlobalStream(RandStream('mt19937ar','seed', 2013));
vv = soundSpeed - soundSpeedWater;
soundSpeed = soundSpeedWater + awgn(vv, homogenityLevel, 'measured');

%% Smooth the phantom with an averaging filter
% reshapedV = reshape(v(:,1), nbTilesX, nbTilesY);
% reshapedV = filter2(averageFilter, reshapedV);
% reshapedV(1:size(averageFilter,1),:) = 1;
% reshapedV(:,1:size(averageFilter,1)) = 1;
% reshapedV(:,end-size(averageFilter,1):end) = 1;
% reshapedV(end-size(averageFilter,1):end,:) = 1;
% v(:,1) = reshapedV(:);

%% Remove noise from water cells.
soundSpeed(waterIndex) = soundSpeedWater;
soundSpeed(ribsIndex) = soundSpeedRibs;

% %%Create the initial sound speed distribution for a given rib phantom.
% soundSpeedInitial = soundSpeed;
% soundSpeedInitial(noRibsIndex) = soundSpeedWater;

% % Build output filename.
% outputDir = [];
% if withRibs == 1, withRibsString = 'withRibs_'; else withRibsString = 'noRibs_'; end
% gridSizeString = [num2str(gridSizeX) '_' num2str(gridSizeY) '_' num2str(gridSizeZ)];
% fileName = [outputDir 'soundSpeedPhantom9_' withRibsString gridSizeString];
% matFileName = [fileName '.mat'];
% matFileNameInitial = [fileName '_initial.mat'];
% binFileName = [fileName '.bin'];
% binFileNameInitial = [fileName '_initial.bin'];
% 
% % Save phantom to a mat file.
% save(matFileName, 'soundSpeed');
% 
% % Save phantom to a binary file.
% save3DArrayToBinaryFile(soundSpeed, binFileName);
% 
% if withRibs == 1
%     % Save the initial distribution for the inversion.
%     save(matFileNameInitial, 'soundSpeedInitial');
%     save3DArrayToBinaryFile(soundSpeedInitial, binFileNameInitial);
% end

%% ----------- END OF CODE --------------
