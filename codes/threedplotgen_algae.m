%clear

%Load Files

baseFileName = 'Mike_data_';
fileExtension = '.dat';
startIndex = 1;
endIndex = 100001;
timesteps = 150;
M = zeros(timesteps,endIndex - startIndex + 1);
for k = 1:timesteps
    fileName = sprintf('%s%02d%s', baseFileName, k, fileExtension);
    B = load(fileName);
    M(k,:) = B(:,2)';
end

baseFileName = 'Ana_data_';
fileExtension = '.dat';
startIndex = 1;
endIndex = 100001;
timesteps = 150;
A = zeros(timesteps,endIndex - startIndex + 1);
for k = 1:timesteps
    fileName = sprintf('%s%02d%s', baseFileName, k, fileExtension);
    B = load(fileName);
    A(k,:) = B(:,2)';
end

baseFileName = 'Cyanototal_data_';
fileExtension = '.dat';
startIndex = 1;
endIndex = 100001;
timesteps = 150;
CT = zeros(timesteps,endIndex - startIndex + 1);
for k = 1:timesteps
    fileName = sprintf('%s%02d%s', baseFileName, k, fileExtension);
    B = load(fileName);
    CT(k,:) = B(:,2)';
end

baseFileName = 'Nitrogen_data_';
fileExtension = '.dat';
startIndex = 1;
endIndex = 100001;
timesteps = 150;
NT = zeros(timesteps,endIndex - startIndex + 1);
for k = 1:timesteps
    fileName = sprintf('%s%02d%s', baseFileName, k, fileExtension);
    B = load(fileName);
    NT(k,:) = B(:,2)';
end

baseFileName = 'Phosphorus_data_';
fileExtension = '.dat';
startIndex = 1;
endIndex = 100001;
timesteps = 150;
PT = zeros(timesteps,endIndex - startIndex + 1);
for k = 1:timesteps
    fileName = sprintf('%s%02d%s', baseFileName, k, fileExtension);
    B = load(fileName);
    PT(k,:) = B(:,2)';
end

baseFileName = 'Chemtotal_data_';
fileExtension = '.dat';
startIndex = 1;
endIndex = 100001;
timesteps = 150;
ChT = zeros(timesteps,endIndex - startIndex + 1);
for k = 1:timesteps
    fileName = sprintf('%s%02d%s', baseFileName, k, fileExtension);
    B = load(fileName);
    ChT(k,:) = B(:,2)';
 end

nt = 150;
xmin = 0;
L = 4e10*10/1e9;  %Domain Length
Ln = 1e5;  %Number of discretized points on the domain
dx = L/Ln;
tmax = timesteps;
[tGrid, xGrid] = meshgrid(linspace(0, tmax, nt), linspace(xmin, L, Ln+1));


figure(1);clf();
h1 = surf(tGrid, xGrid, M', 'EdgeColor','none');
hold on;
%h2 = scatter3(tData, PData, WData, 'ko','filled');
xlabel('t (days)'); ylabel('x (km)'); zlabel('u(x,t) (\mug/L)');
title('Mike Concentration in Time');
% Invert the X and Y axes
set(gca, 'XDir', 'reverse');
% Get the current Axes object
ax = gca;
% Set the FontSize property of the Axes object
ax.FontSize = 20; % You can change 14 to your desired font size
view(45,30); grid on; 
cb = colorbar; % Create colorbar
cb.Location = 'northoutside'; %
exportgraphics(gcf, 'threed_images/Mike_3d_timeplot.png', 'Resolution', 300);
hold off;
drawnow


figure(2);clf();
h1 = surf(tGrid, xGrid, A', 'EdgeColor','none');
hold on;
%h2 = scatter3(tData, PData, WData, 'ko','filled');
xlabel('t (days)'); ylabel('x (km)'); zlabel('w(x,t) (\mug/L)');
title('Ana Concentration in Time');
% Invert the X and Y axes
set(gca, 'XDir', 'reverse');
% Get the current Axes object
ax = gca;
% Set the FontSize property of the Axes object
ax.FontSize = 20; % You can change 14 to your desired font size
view(45,30); grid on; 
cb = colorbar; % Create colorbar
cb.Location = 'northoutside'; %
exportgraphics(gcf, 'threed_images/Ana_3d_timeplot.png', 'Resolution', 300);
hold off;
drawnow


figure(3);clf();
h1 = surf(tGrid, xGrid, CT', 'EdgeColor','none');
hold on;
%h2 = scatter3(tData, PData, WData, 'ko','filled');
xlabel('t (days)'); ylabel('x (km)'); zlabel('u(x,t) + w(x,t) (\mug/L)');
title('Total Cyanobacterial Concentration in Time');
% Invert the X and Y axes
set(gca, 'XDir', 'reverse');
% Get the current Axes object
ax = gca;
% Set the FontSize property of the Axes object
ax.FontSize = 20; % You can change 14 to your desired font size
view(45,30); grid on; 
cb = colorbar; % Create colorbar
cb.Location = 'northoutside'; %
exportgraphics(gcf, 'threed_images/TotalCyano_3d_timeplot.png', 'Resolution', 300);
hold off;
drawnow


figure(4);clf();
h1 = surf(tGrid, xGrid, NT', 'EdgeColor','none');
hold on;
%h2 = scatter3(tData, PData, WData, 'ko','filled');
xlabel('t (days)'); ylabel('x (km)'); zlabel('N(x,t) (\mug/L)');
title('Nitrogen Concentration in Time');
% Invert the X and Y axes
set(gca, 'XDir', 'reverse');
% Get the current Axes object
ax = gca;
% Set the FontSize property of the Axes object
ax.FontSize = 20; % You can change 14 to your desired font size
view(45,30); grid on; 
cb = colorbar; % Create colorbar
cb.Location = 'northoutside'; %
exportgraphics(gcf, 'threed_images/Nitrogen_3d_timeplot.png', 'Resolution', 300);
hold off;
drawnow


figure(5);clf();
h1 = surf(tGrid, xGrid, PT', 'EdgeColor','none');
hold on;
%h2 = scatter3(tData, PData, WData, 'ko','filled');
xlabel('t (days)'); ylabel('x (km)'); zlabel('P(x,t) (\mug/L)');
title('Phosphorus Concentration in Time');
% Invert the X and Y axes
set(gca, 'XDir', 'reverse');
% Get the current Axes object
ax = gca;
% Set the FontSize property of the Axes object
ax.FontSize = 20; % You can change 14 to your desired font size
view(45,30); grid on; 
cb = colorbar; % Create colorbar
cb.Location = 'northoutside'; %
exportgraphics(gcf, 'threed_images/Phosphorus_3d_timeplot.png', 'Resolution', 300);
hold off;
drawnow


figure(6);clf();
h1 = surf(tGrid, xGrid, ChT', 'EdgeColor','none');
hold on;
%h2 = scatter3(tData, PData, WData, 'ko','filled');
xlabel('t (days)'); ylabel('x (km)'); zlabel('N(x,t) + P(x,t) (\mug/L)');
title('Total Chemical Concentration in Time');
% Invert the X and Y axes
set(gca, 'XDir', 'reverse');
% Get the current Axes object
ax = gca;
% Set the FontSize property of the Axes object
ax.FontSize = 20; % You can change 14 to your desired font size
view(45,30); grid on;
cb = colorbar; % Create colorbar
cb.Location = 'northoutside'; %
exportgraphics(gcf, 'threed_images/TotalChem_3d_timeplot.png', 'Resolution', 300);
hold off;
drawnow

figure(7)
subplot(1,2,1)
%plot((1:80)*(1.0/.028)/24,A(:,31250)','o-')
plot((1:150),A(:,31250)','o-')

subplot(1,2,2)

%plot((1:80)*(1.0/.028)/24,PT(:,31250)','o-')
plot((1:150),PT(:,31250)','o-')