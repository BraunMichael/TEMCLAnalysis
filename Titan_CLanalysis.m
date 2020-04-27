clear all; close all
centerWavelength = 735; %in nm
numberWavelengthAverage = 5;
pixelSize = 49; %Size in nm
pixelWindow = 5;
CLmapFileName = sprintf('CL Intensity Map - %i nm', centerWavelength);
distanceFileName = sprintf('AllDistances_%inm', centerWavelength);
firstNNdistanceFileName = sprintf('firstNNDistances_%inm', centerWavelength);
resolution = '-r300'; %Resolution of png output files

%[wavelenthFilename, wavelengthFolderpath] = uigetfile('*.txt', 'Choose the wavelength txt file');
wavelengthFolderpath = 'E:\Google Drive\Research TEM 1-24-2019\SiyingCL\';
wavelenthFilename = 'Spectrum_WavelengthInfo.txt';
cd(wavelengthFolderpath)
rawWavelengths = load(wavelenthFilename);

[dataFilename, datafolderpath] = uigetfile('*.txt', 'Choose the data txt file');
cd(datafolderpath)
rawCLdata = load(dataFilename);
CLSI = [];
for i = 1:length(rawWavelengths)
    CLSI = cat(3, CLSI, rawCLdata(((size(rawCLdata,1)/length(rawWavelengths))*(i-1)+1):(size(rawCLdata,1)/length(rawWavelengths))*i,:));
end

[sortedWavelengths, wavelengthIndices] = sort(abs(bsxfun(@minus, rawWavelengths(:,1), centerWavelength)));
image = mean(CLSI(:,:,wavelengthIndices(1:numberWavelengthAverage)),3);

rawCoords = FastPeakFind(image, 100, (fspecial('gaussian', 7, 1)), 5, 1);
xCoords = rawCoords(1:2:end) * pixelSize;
yCoords = rawCoords(2:2:end) * pixelSize;

imagesc([0 pixelSize*size(image, 2)], [0 pixelSize*size(image, 1)], image); hold on
plot(xCoords, yCoords, 'r+')
title(CLmapFileName,'Interpreter','latex')
xlabel('Distance (nm)','interpreter','latex')
ylabel('Distance (nm)','interpreter','latex')
axis image
hold off

%Have user select points they want to remove
questdlg({'Click and drag to draw a rectangle over the unwanted points' 'A dialog will appear to remove additional points if desired'},'Masking Instructions','OK','Cancel','OK');
rect = getrect;
indicesToRemove = inpolygon(xCoords, yCoords, [rect(1) rect(1) rect(1)+rect(3) rect(1)+rect(3)], [rect(2) rect(2)+rect(4) rect(2)+rect(4) rect(2)]);
xCoords(indicesToRemove) = [];
yCoords(indicesToRemove) = [];

%Select multiple points, using s as break point
s=0;
while s==0
    imagesc([0 pixelSize*size(image, 2)], [0 pixelSize*size(image, 1)], image); hold on
    plot(xCoords, yCoords, 'r+')
    title(CLmapFileName,'Interpreter','latex')
    xlabel('Distance (nm)','interpreter','latex')
    ylabel('Distance (nm)','interpreter','latex')
    axis image
    hold off
    button=questdlg('Do you want to remove more points?','Continue selecting?','yes','no','yes');
    switch button
        case 'yes'
            %Have user select more points they want to remove
            rect = getrect;
            indicesToRemove = inpolygon(xCoords, yCoords, [rect(1) rect(1) rect(1)+rect(3) rect(1)+rect(3)], [rect(2) rect(2)+rect(4) rect(2)+rect(4) rect(2)]);
            xCoords(indicesToRemove) = [];
            yCoords(indicesToRemove) = [];
        case 'no'
            s=1;
    end
end


%% Still working on this
xc = zeros(length(xCoords),1);
yc = zeros(length(xCoords),1);

parfor ii = 1:length(xCoords)
    fprintf('Working on coordinate %i\n', ii);
    [xc(ii), yc(ii)] = gaussian2DImageFit(image, pixelSize, pixelWindow, xCoords(ii), yCoords(ii));
end


% imagesc([0 pixelSize*size(image, 2)], [0 pixelSize*size(image, 1)], image); hold on
% plot(xCoords, yCoords, 'r+')
% plot(xc, yc, 'ko')
% title(CLmapFileName,'Interpreter','latex')
% xlabel('Distance (nm)','interpreter','latex')
% ylabel('Distance (nm)','interpreter','latex')
% axis image
% hold off

%%

if exist(strcat(pwd, '\', CLmapFileName,'.png'), 'file') == 2
    delete(strcat(CLmapFileName,'.png'))
end
print(strcat(CLmapFileName,'.png'),'-dpng',resolution)


%% Calculate Distances
distances = [];
NNdistances = [];
rawNNdistances = [];

for i = 1:length(xCoords)
    for j = 1:length(yCoords)
        if i ~= j
            rawNNdistances = [rawNNdistances sqrt((xCoords(i)-xCoords(j))^2 + (yCoords(i)-yCoords(j))^2)];
        end
    end
    NNdistances = [NNdistances; sort(rawNNdistances)];
    rawNNdistances = [];
end
sortedNNdistances = sort(NNdistances);
for i = 1:length(xCoords)
    if i < length(yCoords)
        for j = (i+1):length(yCoords)
            distances = [distances sqrt((xCoords(i)-xCoords(j))^2 + (yCoords(i)-yCoords(j))^2)];
        end
    end
end

figure()
histogram(distances,20)
title(sprintf('Distances Between Regions - %i nm', centerWavelength), 'Interpreter', 'latex')
xlabel('Distance (nm)','interpreter','latex')
ylabel('Count','interpreter','latex')
axis square

if exist(strcat(pwd, '\', distanceFileName,'.png'), 'file') == 2
    delete(strcat(distanceFileName,'.png'))
end
print(strcat(distanceFileName,'.png'),'-dpng',resolution)


figure()
histogram(NNdistances(:,1),20)
title(sprintf('$$1^{st}$$ Nearest Neighbor Distances - %i nm', centerWavelength),'Interpreter','latex')
xlabel('Distance (nm)','interpreter','latex')
ylabel('Count','interpreter','latex')
axis square

if exist(strcat(pwd, '\', firstNNdistanceFileName,'.png'), 'file') == 2
    delete(strcat(firstNNdistanceFileName,'.png'))
end
print(strcat(firstNNdistanceFileName,'.png'),'-dpng',resolution)

%figure()
% histogram(unique(NNdistances(:,1)),10)
% title('Unique $$1^{st}$$ Nearest Neighbor Distances','Interpreter','latex')
% xlabel('Distance (nm)','interpreter','latex')
% ylabel('Count','interpreter','latex')
% axis square


%% Save data to files
distanceFile = fopen(strcat(distanceFileName, '.txt'),'wt');
fprintf(distanceFile, [ 'Distance (nm)' '\n']);
fprintf(distanceFile, '%f\n', distances');
fclose(distanceFile);

NNdistanceFile = fopen(sprintf('AllNNDistances_%inm.txt', centerWavelength),'wt');
for k = 1:size(NNdistances,2)
    fprintf(NNdistanceFile, '%i NN Distance (nm)\t', k);
end
fprintf(NNdistanceFile, '\n');
for m = 1:size(NNdistances,1)
    fprintf(NNdistanceFile, '%f\t', NNdistances(m,:));
    fprintf(NNdistanceFile, '\n');
end
fclose(NNdistanceFile);

coordFile = fopen(sprintf('Coordinates_%inm.txt', centerWavelength),'wt');
fprintf(coordFile, 'x (pixel)\ty (pixel)\tx (nm)\ty (nm)\n');
fprintf(coordFile, '%i\t%i\n', xCoords./pixelSize, yCoords./pixelSize, xCoords, yCoords);
fclose(coordFile);

% uniqueNNdistanceFile = fopen('UniqueNNDistances.txt','wt');
% fprintf(uniqueNNdistanceFile, [ 'Distance (nm)' '\n']);
% fprintf(uniqueNNdistanceFile, '%f\n', unique(distances)');
