function [distances] = RFP_generator(phi, contrast, amplitude)

% note - this would run faster if saved as a function that you call with RF
% properties and outputs 'distances' variable

%--------------------------------------------------------------------------
% create RF pattern

% set up pix per degree (*think* these values are broadly accurate)
display.width = 52.4;       % Monitor width (cm).
display.horRes = 1920;      % Horizontal resolution of monitor (pixels).
display.viewDist = 57;      % Participant viewing distance (cm).
display.ppd = display.horRes/(2*atand((display.width/2)/display.viewDist)); % Calculates pixels/�. 

%set up RF pattern
RF.numBoundPts = 2^12; % number of boundary points for RF pattern
RF.r0 = 2 * display.ppd; % average radius in degrees vis. angle
% % RF.amp = 0.5;
RF.amp = amplitude;
RF.freq = 3;
% RF.phase = -90 * (pi/180); % phase in deg, convered to rad.
RF.phase = phi; %provide in radians. 
RF.theta = linspace(0,2*pi,RF.numBoundPts + 1);
RF.theta(end) = []; % so no duplicate points
RF.rho = RF.r0*(1 + RF.amp*sin(RF.freq*RF.theta + RF.phase));

% throw RF pattern into boundary with x,y coords
bound = zeros(RF.numBoundPts,2);
[bound(:,1),bound(:,2)] = pol2cart(RF.theta,RF.rho);

%--------------------------------------------------------------------------
% Set up image

% play around with imsize to make sure it's big enough to fully contain
% your largest RF pattern
imSize = 500;

% calculate centroid based on wikipedia formula for centroid of polygon
imArea = (bound(1:end-1,1).*bound(2:end,2))-(bound(2:end,1).*bound(1:end-1,2));
imCent = sum(bsxfun(@times,bound(1:end-1,:)+bound(2:end,:),imArea));
imCent = imCent/(3*sum(imArea));

% centre bound's centroid on image centre
bound = bsxfun(@minus,bound,imCent) + imSize/2;

%--------------------------------------------------------------------------
% Define some key params

% define contrast to be used
% % contrast = 0.5;

% define scaling factor
% Note: keeping this simple and setting scaling factor manually - min it
% should be set to is 1, increases will improve the look of the render at
% the cost of time - for offline generation I typically set it to 10
scaleF = 2;

%{
for a potentially better way to calculate scaling factor - basically take
the average distance between two points on bound (e.g. 0.2px) and use that
to work out max resolution where more or less get adjacent pixels (e.g.
1/0.2 = 5)...add 1 just to up it a bit...
scaleF = ceil(1/mean(hypot(bound(1:end-1,1)-bound(2:end,1),bound(1:end-1,2)-bound(2:end,2))))+1;
%}

% define line thickness (or more specifically, spatial frequency)
lineThickness = 7.5*scaleF;

%--------------------------------------------------------------------------
% Scale up bound and lay it on image

% turn bound (set of x,y co-ordinates) into boundary image
image = false(imSize*scaleF,imSize*scaleF);
bound = round(bound.*scaleF);

% can't remember how this bit works, but basically just puts a 1 wherever
% the RF boundary falls on the image
image(bound(:,2) + (bound(:,1)-1)*size(image,1)) = 1;

%--------------------------------------------------------------------------
% Do the distance transform

% calculate distance transform
distances = bwdist(image,'euclidean');

% resize back down
distances = imresize(distances,1/scaleF);

% adjust to account for eventual line thickness/spatial frequency
distances = distances / lineThickness;

% apply gaussian filter (D4 bit)
distances = (1-4*(distances.^2) + 4/3*(distances.^4)) .* exp(-(distances.^2));

% set contrast
distances = distances .* contrast;

% change allowable range of distances from [-1 +1] to [0 255]
distances = (distances.*127.5)+127.5;

% cast to uint8
distances = uint8(distances);

% rotate image 
distances = flipud(distances);

% % % display the result
% figure; imshow(distances);

