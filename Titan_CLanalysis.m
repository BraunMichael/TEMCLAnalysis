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
wavelenthFilename = 'Spectrum_Wavelengtclear all; close all; warning('off','all');
%set(0,'defaultaxesfontname','arial')
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',32)

pixelSize = 49; %Size in nm
%numberWavelengthAverage = 5;
%pixelWindow = 5; Not averaging at the moment
resolution = '-r300'; %Resolution of png output files
maskFilename = 'spatialMask.mat';
%histogramEdges = (0:25:1000);
pdfX_Values = (0:1:1000);



%[wavelenthFilename, wavelengthFolderpath] = uigetfile('*.txt', 'Choose the wavelength txt file');
wavelengthFolderpath = 'E:\Google Drive\Research TEM\SiyingCL\';
wavelenthFilename = 'Spectrum_WavelengthInfo.txt';
cd(wavelengthFolderpath)
rawWavelengths = load(wavelenthFilename);
wavelengths = rawWavelengths(:,1);




%Set this folder as current directory
datafolderpath=uigetdir('','Select folder with files');
cd(datafolderpath)


tic

D = dir('*.txt');
[tmp, ind] = natsortfiles({D.name});
theFiles = D(ind);

outputFolderpath = strcat(datafolderpath, '\Output');
if exist(outputFolderpath,'dir') ~= 7
    mkdir(outputFolderpath)
end


button=questdlg('Do you want to make a new spatial mask for sample extent? This will delete and update the existing mask','Select new mask?','yes','no','no');
switch button
    case 'yes'
        newmask = 1;
        if exist(strcat(datafolderpath,'\',maskFilename),'file') == 2
            delete(strcat(datafolderpath,'\',maskFilename))
        end
    case 'no'
        newmask = 0;
        load(maskFilename);
end

timeList = [];

for mm = 1:numel(theFiles)
    cd(datafolderpath)
    rawCLdata = load(theFiles(mm).name);
    timeFileCell = regexp(theFiles(mm).name, '\d*','Match');
    timeFile = timeFileCell{1};
    timeList = [timeList str2double(timeFile)];
    CLSI = [];
    for i = 1:length(rawWavelengths)
        CLSI = cat(3, CLSI, rawCLdata(((size(rawCLdata,1)/length(rawWavelengths))*(i-1)+1):(size(rawCLdata,1)/length(rawWavelengths))*i,:));
    end
    image = sum(CLSI,3);
    %centerWavelength = rawWavelengths(maxWavelengthIndex);
    %CLmapFileName = sprintf('CL Intensity Map_%smin_%.0fnm', timeFile, centerWavelength);
    
    %Plot image, select area of interest to create exclusion mask
    figure()
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1);
    imagesc([0 pixelSize*size(image, 2)], [0 pixelSize*size(image, 1)], image);
    title(sprintf('CL Intensity Map_%smin', timeFile),'Interpreter','latex')
    xlabel('Distance (nm)','interpreter','latex')
    ylabel('Distance (nm)','interpreter','latex')
    axis image
    
    
    if newmask
        %Add a ui pop up prompting user to draw exclusion and double click to
        %finish selection
        imageMask = roipoly;
        maskCell{mm} = imageMask;
        save(maskFilename,'maskCell')
    else
        imageMask = maskCell{mm};
    end
    
    
    cd(outputFolderpath)
    [~, maxImageIndicies] = max(CLSI,[],3);
    maxWavelengthImage = wavelengths(maxImageIndicies);
    
    %wavelengthToComposition = 0.5867*lambda - 310.725 %Linear fit of data
    %provided by Siying. Gives Iodine %
    compositionImage = (0.5867.*maxWavelengthImage - 310.725) .* imageMask;
    figure()
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1);
    imagesc([0 pixelSize*size(image, 2)], [0 pixelSize*size(image, 1)], compositionImage);
    title(sprintf('Composition Map_%smin', timeFile),'Interpreter','latex')
    xlabel('Distance (nm)','interpreter','latex')
    ylabel('Distance (nm)','interpreter','latex')
    axis image
    colormap(jet(4096))
    colorbar
    pause(0.5);
    
    CompMapFileName = sprintf('Composition Map_%smin', timeFile);
    if exist(strcat(pwd, '\', CompMapFileName,'.png'), 'file') == 2
        delete(strcat(CompMapFileName,'.png'))
    end
    print(strcat(CompMapFileName,'.png'),'-dpng',resolution)
    integratedIntensity = permute(sum(sum(bsxfun(@times, CLSI, cast(imageMask, class(CLSI))),1),2),[3 2 1]);
    normalizedIntegratedIntensity = integratedIntensity./max(integratedIntensity);
    figure()
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1);
    plot(wavelengths, normalizedIntegratedIntensity, 'k')
    title(sprintf('Normalized Integrated Spectrum_%smin', timeFile),'Interpreter','latex')
    xlabel('Wavelength (nm)','interpreter','latex')
    ylabel('Normalized Integrated Intensity','interpreter','latex')
    pause(0.5);
    
    integratedSpectraFileName = sprintf('IntegratedSpectrum_%smin', timeFile);
    if exist(strcat(pwd, '\', integratedSpectraFileName,'.png'), 'file') == 2
        delete(strcat(integratedSpectraFileName,'.png'))
    end
    print(strcat(integratedSpectraFileName,'.png'),'-dpng',resolution)
    integratedSpectraDataFile = fopen(sprintf('IntegratedSpectra_%smin.txt', timeFile),'wt');
    fprintf(integratedSpectraDataFile, 'Wavelength (nm)\tIntegrated Intensity (arb units)\tNormalized Integrated Intensity (arb units)\n');
    for ii = 1:length(wavelengths)
        fprintf(integratedSpectraDataFile, '%f\t%f\t%.8f\n', wavelengths(ii), integratedIntensity(ii), normalizedIntegratedIntensity(ii));
    end
    fclose(integratedSpectraDataFile);
    
    
    tic
    parfor nn = 1:length(rawWavelengths)
        %Make all values outside selected area 0
        image = CLSI(:,:,nn) .* imageMask;
        centerWavelength = rawWavelengths(nn);
        CLmapFileName = sprintf('CL Intensity Map_%smin_%.0fnm', timeFile, centerWavelength);
        NNdistanceFileName = sprintf('NNDistances_%smin_%.0fnm', timeFile, centerWavelength);
        distanceFileName = sprintf('AllDistances_%smin_%.0fnm', timeFile, centerWavelength);
        
        
        
        
        
        rawCoords = FastPeakFind(image, 100, (fspecial('gaussian', 7, 1)), 5, 1);
        xCoords = rawCoords(1:2:end) * pixelSize;
        yCoords = rawCoords(2:2:end) * pixelSize;
        
        figure()
        pause(0.00001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1);
        imagesc([0 pixelSize*size(image, 2)], [0 pixelSize*size(image, 1)], image); hold on
        plot(xCoords, yCoords, 'r+')
        title(CLmapFileName,'Interpreter','latex')
        xlabel('Distance (nm)','interpreter','latex')
        ylabel('Distance (nm)','interpreter','latex')
        axis image
        hold off
        
        
        %% Still working on this
        % xc = zeros(length(xCoords),1);
        % yc = zeros(length(xCoords),1);
        %
        % parfor ii = 1:length(xCoords)
        %     fprintf('Working on coordinate %.0f\n', ii);
        %     [xc(ii), yc(ii)] = gaussian2DImageFit(image, pixelSize, pixelWindow, xCoords(ii), yCoords(ii));
        % end
        
        
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
        
        if ~isempty(xCoords)
            coordFile = fopen(sprintf('Coordinates_%smin_%.0fnm.txt', timeFile, centerWavelength),'wt');
            fprintf(coordFile, 'x (pixel)\ty (pixel)\tx (nm)\ty (nm)\n');
            for qq = 1:length(xCoords)
                fprintf(coordFile, '%.0f\t%.0f\t%.0f\t%.0f\n', xCoords(qq)./pixelSize, yCoords(qq)./pixelSize, xCoords(qq), yCoords(qq));
            end
            fclose(coordFile);
            if length(xCoords) > 1
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
                
                
                %% Save data to files
                distanceFile = fopen(strcat(distanceFileName, '.txt'),'wt');
                fprintf(distanceFile, [ 'Distance (nm)' '\n']);
                fprintf(distanceFile, '%f\n', distances');
                fclose(distanceFile);
                
                NNdistanceFile = fopen(strcat(NNdistanceFileName, '.txt'),'wt');
                for k = 1:size(NNdistances,2)
                    fprintf(NNdistanceFile, '%.0f NN Distance (nm)\t', k);
                end
                fprintf(NNdistanceFile, '\n');
                for m = 1:size(NNdistances,1)
                    fprintf(NNdistanceFile, '%f\t', NNdistances(m,:));
                    fprintf(NNdistanceFile, '\n');
                end
                fclose(NNdistanceFile);
            end
        end
        if length(xCoords) > 15
            figure()
            hold on
            probDistribution = fitdist(NNdistances(:,1),'Lognormal');
            mu = probDistribution.mu;
            sigma = probDistribution.sigma;
            pdfY_Values = pdf(probDistribution, pdfX_Values);
            probDistributionNormal = fitdist(NNdistances(:,1),'Normal');
            title(strcat('First ', NNdistanceFileName),'Interpreter','latex')
            plot(pdfX_Values, pdfY_Values, 'r', pdfX_Values, pdfY_Values,'b')
            xlabel('Distance (nm)','interpreter','latex')
            ylabel('Probability','interpreter','latex')
            xlim([0 1000])
            yLimits = ylim;
            ylim([0 yLimits(2)])
            axis square
            stemY = 0.1.*yLimits(2).*ones(length(NNdistances(:,1)),1);
            stem(NNdistances(:,1), stemY,'Marker','none','LineWidth',1,'Color','black')
            hold off
            
            if exist(strcat(pwd, '\first', NNdistanceFileName,'.png'), 'file') == 2
                delete(strcat('first', NNdistanceFileName,'.png'))
            end
            print(strcat('first', NNdistanceFileName,'.png'),'-dpng',resolution)
            
            [~, modeIndex] = max(pdf(probDistribution, pdfX_Values));
            modeMatrix(mm,nn) = pdfX_Values(modeIndex);
        else
            modeMatrix(mm,nn) = 0;
            mu = 0;
            sigma = 0;
            pdfY_Values = zeros(size(pdfX_Values));
        end
        
        pdfMuList(nn,mm) = mu;
        pdfSigmaList(nn,mm) = sigma;
        pdfYValueList(:,nn,mm) = pdfY_Values;
        
    end
    toc
    tempWavelengths = [];
    tempPdfYValueList = [];
    for ii = 1:length(wavelengths)
        if any(pdfYValueList(:,ii,mm))
            tempPdfYValueList = [tempPdfYValueList pdfYValueList(:,ii,mm)];
            tempWavelengths = [tempWavelengths; wavelengths(ii)];
        end
    end
    
    
    pdfYValueDataFile_headerstring = 'Distance (nm)\t';
    for ii = 1:(length(tempWavelengths)-1)
            pdfYValueDataFile_headerstring = strcat(pdfYValueDataFile_headerstring,sprintf('Probability at %.0fnm\\t',tempWavelengths(ii)));
    end
    pdfYValueDataFile_headerstring = strcat(pdfYValueDataFile_headerstring,sprintf('Probability at %.0fnm\\n',tempWavelengths(ii+1)));
    
    
    pdfYValueDataFile = fopen(sprintf('pdfYValueDataFile_%smin.txt', timeFile),'wt');
    fprintf(pdfYValueDataFile, pdfYValueDataFile_headerstring);
    for ii = 1:length(pdfX_Values)
        fprintf(pdfYValueDataFile, '%.0f\t', pdfX_Values(ii));
        fprintf(pdfYValueDataFile, '%e\t', tempPdfYValueList(ii,:));
        fprintf(pdfYValueDataFile, '\n');
    end
    fclose(pdfYValueDataFile);
    
end




pdfParamDataFile_headerstring = 'Wavelength (nm)\t';
for ii = 1:(length(timeList)-1)
    pdfParamDataFile_headerstring = strcat(pdfParamDataFile_headerstring,sprintf('mu at %.0fmin\\tsigma at %.0fmin\\t',timeList(ii),timeList(ii)));
end
pdfParamDataFile_headerstring = strcat(pdfParamDataFile_headerstring,sprintf('mu at %.0fmin\\tsigma at %.0fmin\\n',timeList(ii+1),timeList(ii+1)));

for ii = 1:length(timeList)
    pdfCombinedParams(:,2*ii-1) = pdfMuList(:,ii);
    pdfCombinedParams(:,2*ii) = pdfSigmaList(:,ii);
end

pdfParamDataFile = fopen('LognormalProbabilityFunctionParameters.txt','wt');
fprintf(pdfParamDataFile, pdfParamDataFile_headerstring);
for ii = 1:length(wavelengths)
    fprintf(pdfParamDataFile, '%f\t%f\t', wavelengths(ii), pdfCombinedParams(ii,:));
    fprintf(pdfParamDataFile, '\n');
end

fclose(pdfParamDataFile);


cd(datafolderpath)
save('SummaryVariables.mat','modeMatrix','pdfMuList','pdfSigmaList','pdfYValueList','wavelengths')





figure()
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
imagesc([wavelengths(1) wavelengths(end)],[timeList(1) timeList(end)],modeMatrix)
xlabel('Wavelength (nm)','interpreter','latex')
ylabel('Time (min)','interpreter','latex')
h = colorbar;
set(get(h,'label'),'string','First Nearest Neighbor Distance (nm)');
%colormap(b2r(4096))
colormap(hot)hInfo.txt';
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
