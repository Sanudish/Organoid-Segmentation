Debug = 0;
    
[fileA,pathA] = uigetfile('*.tif','SELECT BRIGHTFIELD IMAGE','SELECT BRIGHTFIELD IMAGE');
disp('Brightfield image selected')
[fileAF,pathAF] = uigetfile('*.tif','SELECT FLUORESCENT IMAGE','SELECT FLUORESCENT IMAGE');
disp('Fluorescent image selected')
DESTINATION_FOLDER = uigetdir('C:\','SELECT DESTINATION FOLDER');
disp('Destination folder selected')

prompt = {'Enter Mean Intensity Threshold:','Enter Area Threshold:','Enter Background Threshold:','Enter Edge Sensitivity'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'40','20','15','0.5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

Thresh = str2num(answer{1});
AThresh = str2num(answer{2});
BThresh = str2num(answer{3});
Sensitivity = str2num(answer{4});

A = imread(fullfile(pathA,fileA));
AF = imread(fullfile(pathAF,fileAF));

disp('Closing open figures')
close all

if Debug == 1
    disp('Generating brightfield image')
    %Plot original bright field image at 0hr
    figure
    imshow(A),colorbar,caxis([0,3000]),title('Brightfield Image')
end 


%Identify threshold for sobel and add adjusting factor
[~,thresholdA] = edge(A, 'sobel');
fudgeFactor = Sensitivity;
BWsA = edge(A,'sobel',thresholdA*fudgeFactor);

se90A = strel('line',3,90);
se0A = strel('line',3,0);

BWsdilA = imdilate(BWsA,[se90A,se0A]);
BWdfillA = imfill(BWsdilA, 'holes');

seDA = strel('diamond',1);
BWfinalA = imerode(BWdfillA, seDA);
BWfinalA = imerode(BWfinalA, seDA);
disp('Finished Sobel Segmentation')
%%%%%%%%%%%%%%%%% Preparing sample for watershed algorithm
Sobelsegmentation = BWfinalA;
segmentopen = ~bwareaopen(~Sobelsegmentation,10);

%%%
D = -bwdist(~Sobelsegmentation);
%%%
disp('First Watershed')
FirstWS = watershed(D);
%%%

segmentopen = Sobelsegmentation;
segmentopen(FirstWS == 0) = 0;
%%%

mask = imextendedmin(D,2);
%%%
disp('Second Watershed')
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = Sobelsegmentation;
bw3(Ld2 == 0) = 0;
%%%%%%%%%%%%%%%%%%%

if Debug == 1
    disp('Generating segmented brightfield image')
    % Plot the segmented brightfield Image at 0hr
    figure
    imshow(bw3),colorbar,caxis([0,1]),title('Segmented Brightfield Image')

    disp('Generating fluorescent image')
    %Plot the fluorscent image at 0hr
    figure
    imshow(AF),colorbar,caxis([0,800]),title('Fluorescent Image')
    
    disp('Generating fluorescent image - with mask')
    %Highlight the regions at which we will be able to see through by setting
    %the segmented region to max intensity on scale 
    burnedimage = AF;
    burnedimage(bw3) = 800;
    figure
    imshow(burnedimage),colorbar,caxis([0,800]),title('Fluorescent Image - With Mask')

    disp('Generating masked image - black outside mask, fluorescent inside')
    %Apply mask over fluorscent image to effectively see through the segmented
    %regions and plot
    blackMaskedImage = AF;
    blackMaskedImage(~bw3) = 0;
    figure
    imshow(blackMaskedImage),colorbar,caxis([0,800]),title('Masked Image - Black Outside Mask, Fluorescent Inside');
    
    disp('Generating masked image with fluorescent spots inside')
    %Apply mask over fluorscent image to effectively see through the segmented
    %regions and plot
    SpotsInMask = AF;
    SpotsInMask(~bw3) = 10000;
    figure
    imshow(SpotsInMask),colormap(turbo(256)),colorbar,caxis([0,800]),title('Masked Image - Fluorescence Inside');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now use bwconncomp function to identify connected components in masked
%image - because mask is applied, the area around high intensity pixels are
%considered

%Find area and meanintensity, using labelmatrix to index masked image
blackMaskedImage = AF;
blackMaskedImage(~bw3) = 0;
CC = bwconncomp(blackMaskedImage);
L = labelmatrix(CC);
stats = regionprops('table',L,blackMaskedImage,'Area','MeanIntensity');

disp('Generating Table with Area and Mean Intensities')
%Generate table with above desired threshold for mean intensity
numstats = table2array(stats);
label = find(numstats(:,2)>Thresh);
f = transpose(1:length(label));
AMI = numstats(label,1:2);
combine = [f,AMI];

%Filter for area threshold
labelarea = find(combine(:,2)>AThresh);
g = transpose(1:length(labelarea));
AAMI = combine(labelarea,2:3);
Finalcombine = [g,AAMI];

BorF = Finalcombine(:,3);
BorF(BorF>=BThresh)=70;
BorF(BorF<BThresh) = 66;
BorFletter = char(BorF);
BorFstr = cellstr(BorFletter);
FINALTABLE = table(Finalcombine(:,1),Finalcombine(:,2),Finalcombine(:,3),BorFstr,'VariableNames',{'Label','Area(Px)','MeanIntensity','B or F'});

writetable(FINALTABLE,fullfile(DESTINATION_FOLDER,'Area&MeanIntensity.csv'));


%Generate figure for areas which are plotted

disp('Generating areas for which stats are calculated')
areas1 = ismember(L,label);

%Convert logical label matrix to sequential numbers
Sequence = bwlabel(areas1);

areas2 = ismember(Sequence,labelarea);

if Debug == 1
    figure
    imshow(areas2),colorbar,title('Areas for which stats are calculated.');
end

disp('Generating Labels')
%Generate labels
sequence2 = bwlabel(areas2);
finalstats = regionprops('table',sequence2,'Centroid');
finalstatsarr = table2array(finalstats);

x = finalstatsarr(:,1);
y = finalstatsarr(:,2);

if Debug == 1
    figure
    imshow(sequence2)
    hold on
        for ii = 1:length(x)
             text(x(ii),y(ii),num2str(ii),'Color','m','FontSize',14)
        end
end

disp('Generating labelled segmented region over brightfield')
figure
imshow(A),colorbar,caxis([0 3000])
hold on        
        for ii = 1:length(x)
             text(x(ii),y(ii),num2str(ii),'Color','m','FontSize',14)
        end
[Bla,Lla] = bwboundaries(sequence2,'noholes');
for k = 1:length(Bla)
   boundary = Bla{k};
   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1)
end
saveas(figure(1),fullfile(DESTINATION_FOLDER,'1_Labelled_Segmented_Region_Over_Brightfield.png'))

disp('Generating labelled segmented region over fluorescence')
figure
imshow(AF),colorbar,caxis([0 800])
hold on        
        for ii = 1:length(x)
             text(x(ii),y(ii),num2str(ii),'Color','m','FontSize',14)
        end
for k = 1:length(Bla)
   boundary = Bla{k};
   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1)
end
saveas(figure(2),fullfile(DESTINATION_FOLDER,'2_Labelled_Segmented_Region_Over_Fluorescence.png'))


disp('end')



    



