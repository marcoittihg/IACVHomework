%%LOAD IMAGE
im = imread('First.jpg');
baseImage = double(im)/255;
bwImage = rgb2gray(baseImage);
bwImage = histeq(bwImage);

%% EDGE DETECTION
sigma = 2;
sx=[1 0 -1 ; 2 0 -2; 1 0 -1];
sy=sx';

% Convolution with gaussian filter
gSmoothed = imgaussfilt(bwImage, sigma);

% Compute gradient magnitude and direction
dx = conv2(gSmoothed, sx, 'same');
dy = conv2(gSmoothed, sy, 'same');
mag = sqrt(dx.^2+dy.^2);
theta = atand(dy./(dx+0.01));

% NON MAXIMUM SUPPRESSION
sze = size(theta);
outDir = zeros(sze(1), sze(2), 2);
for i = 1:sze(1)
    for j = 1:sze(2)
        dir = theta(i,j);
        v = dir / 90;

        if abs(v) > 3/4
            val = [1 0];
        elseif v >= 1/4 && v <= 3/4
            val = [1 1];
        elseif v >= -3/4 && v <= -1/4
            val = [-1 1];
        else
            val = [0 1];
        end
        
        outDir(i,j,1) = val(1);
        outDir(i,j,2) = val(2);
    end
end

outMag = mag;
for i = 2:sze(1)-1
    for j = 2:sze(2)-1
        
        p = mag(i,j);
        pn = mag(i+outDir(i,j,1), j+outDir(i,j,2));
        pp = mag(i-outDir(i,j,1), j-outDir(i,j,2));
        
        if p < pn || p < pp
            outMag(i,j) = 0;
        end
    end
end

% Thresholding
hTh = 0.1;
lTh = hTh * 0.4;

for i = 1:sze(1)
    for j = 1:sze(2)
        if outMag(i,j) > hTh
            outMag(i,j) = 1;
        elseif outMag(i,j) < lTh
            outMag(i,j) = 0;
        else
            outMag(i,j) = 0.5;
        end
    end
end

found = 1;
while found == 1
    found = 0;
    
    for i = 2:sze(1)-1
        for j = 2:sze(2)-1
            if outMag(i,j) == 0.5
                if outMag(i+1,j) == 1 || outMag(i+1,j+1) == 1 || outMag(i,j+1) == 1 || outMag(i-1,j+1) == 1 || outMag(i-1,j) == 1 || outMag(i-1,j-1) == 1 || outMag(i,j-1) == 1 || outMag(i+1,j-1) == 1
                    outMag(i,j) = 1;
                    found = 1;
                end
            end
        end
    end
end

for i = 2:sze(1)-1
    for j = 2:sze(2)-1
        if outMag(i,j) == 0.5
            outMag(i,j) = 0;
        end
    end
end


%matCanny = edge(bwImage,'Canny');
% imshow([outMag matCanny],[]);
imshow(outMag,[]);

%% CORNERS
kSize = 10;
sigma = 2;
sx=[1 0 -1 ; 2 0 -2; 1 0 -1];
sy=sx';

dx = conv2(bwImage, sx, 'same');
dy = conv2(bwImage, sy, 'same');
dxx = dx.*dx;
dyy = dy.*dy;
dxy = dx.*dy;

% Compute matrix M
sze = size(bwImage);
M = zeros(sze(1), sze(2), 2, 2);
M(:,:,1,1) = imgaussfilt(dxx,2);
M(:,:,2,2) = imgaussfilt(dyy,2);
M(:,:,2,1) = imgaussfilt(dxy,2);
M(:,:,1,2) = imgaussfilt(dxy,2);

% Compute CM
det = M(:,:,1,1).*M(:,:,2,2) - M(:,:,1,2).*M(:,:,2,1);
tr = M(:,:,1,1)+M(:,:,2,2);
CM = det./(tr+0.001);

% Find local maxima
regionalmax = reshape(imregionalmax(CM) .* CM,[],1);
regMax = maxk(regionalmax, 500);

disp(regMax);
regionalmax = imregionalmax(CM);
th = regMax(end);
img = bwImage;
imshow(img);
hold on;
for x = 1:sze(1)
    for y = 1:sze(2)
        if CM(x,y) >= th && regionalmax(x,y) == 1
            disp("Plot");
        end
    end
end
%% Comparison with std implementation
pts = detectHarrisFeatures(bwImage);
plot(pts.selectStrongest(500));

%% LINES
%Data points
count =0;
for x=1:3096
    for y=1:4128
        if outMag(x,y)==0
            continue;
        end
        count=count+1;
    end
end
DPts = zeros(count, 2);
count =0;
for x=1:3096
    for y=1:4128
        if outMag(x,y)==0
            continue;
        end
        DPts(count+1,1) = x;
        DPts(count+1,2) = y;
        
        count=count+1;
    end
end

% Hough Transform
thr=1;
maxRho = norm(size(outMag));
maxTheta = 90;
Model = zeros(2*maxTheta+2, 2*maxRho+2);

dpCount = size(DPts);
dpCount = dpCount(1);

thetaValues = zeros([2*maxTheta]);
for theta=-maxTheta:maxTheta
    thetaValues(theta+maxTheta+1) = theta;
end
cosTheta= cosd(thetaValues);
sinTheta= sind(thetaValues);

% Foreach data point add a vote to all compatible models
for iDp=1:dpCount
    if mod(iDp, 10000) == 0
        disp([iDp dpCount]);
    end
    
    x = DPts(iDp,1);
    y = DPts(iDp,2);
    
    min = -maxRho;
    max = maxRho;
    
    for theta=-maxTheta:maxTheta
        th = theta;
        rho = x*cosTheta(theta+maxTheta+1)+y*sinTheta(theta+maxTheta+1);
        rhoIdx = round(rho);
        
        idx = rhoIdx;
        found = 0;
        while found==0
            if rho-idx < thr
                Model(th+maxTheta+1, idx+maxRho+1) = Model(th+maxTheta+1, idx+maxRho+1)+1;
                idx=idx-1;
            else
                found = 1;
            end
        end
        
        idx = rhoIdx+1;
        found = 0;
        while found==0
            if idx-rho < thr
                Model(th+maxTheta+1, idx+maxRho+1) = Model(th+maxTheta+1, idx+maxRho+1)+1;
                idx=idx+1;
            else
                found = 1;
            end
        end
    end
end

%
maxModels = imregionalmax(Model);
cont=0;
for i=1:181
    for j=1:10321
       cont = cont + maxModels (i,j);
    end
end

hold off
imshow(outMag,[]);
hold on

% Display lines with at least thVotes votes
thVotes = 1500;
lines = cont;
cont=0;
for i=1:181
    for j=1:10321
       if(maxModels (i,j) == 0)
           continue;
       end
       
       cont=cont+1;
       if mod(cont, 100) == 0
           disp([cont lines]);
       end
       
       if Model(i,j) < thVotes
           continue;
       end
       
       theta = i-maxTheta-1;
       rho = j-maxRho-1;
       
       
       x0 = rho/sind(theta);
       y0 = rho/cosd(theta);
       xf = (rho-4128*cosd(theta))/sind(theta);
       yf = (rho-3096*sind(theta))/cosd(theta);
       
       if abs(x0) == Inf
           startL = [y0 0];
           endL = [yf 0];
       elseif abs(y0) == Inf
           startL = [0 x0];
           endL = [0 xf];
       else
           startL = [0 x0];
           endL = [4128 xf];
       end
       
       line([startL(2) endL(2)],[startL(1) endL(1)],'Color','red','LineWidth',2);
    end
end
