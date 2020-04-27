function [xc, yc] = gaussian2DImageFit(image, pixelSize, pixelWindow, xCoord, yCoord)


image_x = (0:1:size(image, 2))*pixelSize;
image_x(end) = [];
image_y = (0:1:size(image, 1))*pixelSize;
image_y(end) = [];
imageY = image_y(((yCoord/pixelSize)-pixelWindow):((yCoord/pixelSize)+pixelWindow));
imageX = image_x(((xCoord/pixelSize)-pixelWindow):((xCoord/pixelSize)+pixelWindow));
S = image(((yCoord/pixelSize)-pixelWindow):((yCoord/pixelSize)+pixelWindow),((xCoord/pixelSize)-pixelWindow):((xCoord/pixelSize)+pixelWindow));


[xData, yData, zData] = prepareSurfaceData(imageX, imageY, S);
xyData = {xData,yData};

[amp, ind] = max(zData); % amp is the amplitude.
xo = xCoord; % guess that it is at the maximum
yo = yCoord; % guess that it is at the maximum
ang = 45; % angle in degrees.
sy = (max(yData)-min(yData))/2;
sx = (max(xData)-min(xData))/2;
zo = median(zData)-1.5*std(zData);
xmax = max(xData)+2;
ymax = max(yData)+2;
xmin = min(xData)-2;
ymin = min(yData)-2;

%% Set up fittype and options.
StartPoint = [amp, ang, sx, sy, xo, yo, zo];%[amp, angle, sx, sy, xo, yo, zo];
Lower = [0.1*amp, 0, 0, 0, xmin, ymin, 0]; %[amp, angle, sx, sy, xo, yo, zo];
Upper = [2*amp, 180, Inf, Inf, xmax, ymax, amp]; % angles greater than 90 are redundant


tols = 1e-16;
%options = optimset('Algorithm','levenberg-marquardt','Display','off','MaxFunEvals',5e3,...
%    'MaxIter',5e3,'TolX',tols,'TolFun',tols,'TolCon',tols ,'UseParallel',false);

options = optimset('Display','off','MaxFunEvals',5e3,...
    'MaxIter',5e3,'TolX',tols,'TolFun',tols,'TolCon',tols ,'UseParallel',false);
%% perform the fitting
gaussian2D = @(par,xy) par(7) + par(1)*exp(-(((xy{1}-par(5)).*cosd(par(2))+(xy{2}-par(6)).*sind(par(2)))./par(3)).^2-((-(xy{1}-par(5)).*sind(par(2))+(xy{2}-par(6)).*cosd(par(2)))./par(4)).^2);
[fitresult,~,~] = lsqcurvefit(gaussian2D,StartPoint,xyData,zData,Lower,Upper,options);
%[fitresult,~,~] = lsqcurvefit(gaussian2D,StartPoint,xyData,zData,[],[],options);

xc = fitresult(5);
yc = fitresult(6);
