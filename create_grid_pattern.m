function [map] = create_grid_pattern(orientation,spacing,phase,boxSize,peakFR,fieldSize)


%% create a grid pattern from sine waves 

% code is adapted from:
% http://www.icn.ucl.ac.uk/courses/MATLAB-Tutorials/Elliot_Freeman/html/gabor_tutorial.html

% due to numerical things, it is best to take double the box size at first,
% and then take the center of the grid pattern
boxSize = boxSize*2; 
imSize = boxSize*2;       % image size: n X n
% boxSize = 50;                         % size of environment - 50 cm x 50 cm
pixelPerCm = imSize/boxSize;            % pixel/cm (spacing is cm/cycle)              
wavelength = spacing*pixelPerCm;        % wavelength (number of pixels/ cycle)
% orientation = 15;                     % grating orientation
% phase = .25;                          % phase (0 -> 1)
X = 1:imSize;                           % X is a vector from 1 to imageSize
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5 
freq = imSize/wavelength;               % compute frequency from wavelength
phaseRad = (phase * 2* pi);             % convert to radians: 0 -> 2*pi

% create the first wave
orientation1 = orientation;
thetaRad1 = (orientation1 / 360) * 2*pi; % convert theta (orientation) to radians
[Xm, Ym] = meshgrid(X0, X0);            % 2D matrices
Xt = Xm * cos(thetaRad1);                % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad1);                % compute proportion of Ym for given orientation
XYt = Xt + Yt;                          % sum X and Y components
XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
grating1 = sin( XYf + phaseRad);        % make 2D sinewave
lines1 = (grating1 > 0.99);

% create the second wave
orientation2 = orientation+60;
thetaRad2 = (orientation2 / 360) * 2*pi;  % convert theta (orientation) to radians
[Xm, Ym] = meshgrid(X0, X0);             % 2D matrices
Xt = Xm * cos(thetaRad2);                % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad2);                % compute proportion of Ym for given orientation
XYt = Xt + Yt;                          % sum X and Y components
XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
phaseRad2 = (rand * 2* pi);  
grating2 = sin( XYf + phaseRad2);        % make 2D sinewave
lines2 = grating2 > 0.99;

% create the third wave
orientation3 = orientation+120;
thetaRad3 = (orientation3 / 360) * 2*pi;  % convert theta (orientation) to radians
[Xm, Ym] = meshgrid(X0, X0);             % 2D matrices
Xt = Xm * cos(thetaRad3);                % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad3);                % compute proportion of Ym for given orientation
XYt = Xt + Yt;                          % sum X and Y components
XYf = XYt * freq * 2*pi; % convert to radians and scale by frequency

% find the phase shift that will give constructive interference for all 3
% sine waves
[~,ind] = max(grating1(:)+grating2(:));
[i,j] = ind2sub(size(Xm),ind);
y = Ym(i,1); x = Xm(1,j);
phaseRad3 = asin(sin((x*cos(thetaRad2)+y*sin(thetaRad2))*freq*2*pi+phaseRad2)) ...
    - (x*cos(thetaRad3)+y*sin(thetaRad3))*freq*2*pi;

grating3 = sin( XYf + phaseRad3);                   % make 2D sinewave
lines3 = grating3 > 0.99;

%{
figure()
subplot(1,4,1)
imagesc( lines1);            
axis off; axis image; 

subplot(1,4,2)
imagesc(lines2);            
axis off; axis image; 

subplot(1,4,3)
imagesc(lines3);            
axis off; axis image; 

subplot(1,4,4)
imagesc( lines1 + lines2 + lines3);            
axis off; axis image;                   % use gray colormap

keyboard
%}

sine_map = grating1 + grating2 + grating3;

%% locate the centers of the fields

% find all the centroids
s = regionprops(sine_map >= 0.99*max(max(sine_map)),'centroid');
centroids = cat(1, s.Centroid);

% plot if I want
%{
subplot(1,3,1)
imagesc(sine_map);
hold on
plot(centroids(:,1),centroids(:,2), 'r+')
hold off
axis off; axis image
%}

%% put a Gaussian blob at each centroid
% use the provided peakFR and the field size

gridPattern = zeros(imSize);
Sigma = fieldSize/(-2*log(0.2)*pi)*eye(2);
[X1,X2] = meshgrid(linspace(0,boxSize,imSize),linspace(0,boxSize,imSize));

centroids = round(centroids);

for k = 1:size(centroids,1)
    mu = [X1(1,centroids(k,1)) X2(centroids(k,2),1)];
    field = reshape(mvnpdf([X1(:) X2(:)],mu,Sigma),imSize,imSize);
    field = peakFR.*field./max(max(field));
    gridPattern = gridPattern+field;
    
end
%{
subplot(1,3,2)
imagesc(gridPattern)
axis off; axis image
%}

% take the inner part of the grid pattern
% inside = imSize/4-imSize/6+1:imSize/2+imSize/6;
inside = imSize/4+1:3*imSize/4;
gridPattern_inside = gridPattern(inside,inside);
map = gridPattern_inside;

%{
subplot(1,3,3)
imagesc(map)
axis off; axis image
%}

return