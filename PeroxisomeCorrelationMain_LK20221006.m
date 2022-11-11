%% PeroxisomeCorrelationMain_LK20221006.m
% image analysis to perform image correlation on mCherry and GFP-labeled
% peroxisome proteins
% Analysis includes sections to:
% 1. Load data (read comments in this section about the input data format)
% 2. Image registration of GFP and mCherry channels based on control point selection,
% correlation, and image transformation using fiducial bead markers
% 3. Find location of peroxisomes using existing particle tracking code
% (requires Troika single particle tracking code maintained by Prof.
% Christy Landes' group available at:
% https://github.com/LandesLab/Troika-Single-particle-tracking; relevant
% citation is DOI: 10.1039/C3CP53968G
% 4. For each peroxisome, perform image fluorscence cross correlation
% spectroscopy at the 17 pixels surrounding the centroid location (read
% comments in this section regarding format of the data saved)
% 5. Save data (optional) 
% 6. Make figures of results (optional)
% Contact Prof. Lydia Kisley, Case Western Reserve University
% lydia.kisley@case.edu for any questions

close all
clear all

%% Peroxisome Analysis %%
% LK20221006
% Remove unused code and comment for manuscript review
% LK20211115
% Condense Zhenghao's "Peroxisome-main" scripts into a single
% script so no need to repeat many selections

% Need to add script here to cycle through all folders and files %

%% 1. Load data
% Data is a .mat 2x4 cell with the following format:
% Variable name is 'Data'
% Row 1: channel/dataset names: ['GFP', 'mCherry', 'bead', 'Cy5']
% Row 2: [1200 x 1200 x 150 double of GFP channel, 1200 x 1200 x 150 double
% of mCherry, 1200 x 2400 x 1 unint16 frame fiducial bead markers, 
% 1200 x 2400 x 1 unint16 frame fiducial bead markers taken with Cy5
% filter] 
% The Cy5 data does not end up being used
fpath='C:\...\'; %folder directory of data
fname='90minData'; %name of .mat file containing 'Data'
saveyn=1; %save the result of not

load(strcat(fpath,fname,'.mat'));
 
%% 2. Open the alignment image with beads
%split into two images
%convert to mat2gray for cpselect
beadleft=mat2gray(Data{2,3}(:,1:size(Data{2,3},2)/2)); %green channel to be moved
beadright=mat2gray(Data{2,3}(:,(size(Data{2,3},2)/2)+1:size(Data{2,3},2))); %red channel to be fixed

%cpselect
[movingPoints,fixedPoints] = cpselect(beadleft,beadright,'Wait',true) 
% cpselect(beadleft,beadright); %need to figure out how to get this to run again after selecting points
movingPointsAdjusted = cpcorr(movingPoints,fixedPoints,beadleft,beadright) %fine-tune locations with correlation
clear beadleft beadright

%apply transformation to GFP and mCherry data ??? 
tform = fitgeotrans(movingPointsAdjusted,fixedPoints,'affine');
sizeim=[size(Data{2,1},1),size(Data{2,1},2)]; %get size, so transform doesn't change  matrix size
framesg_reg = imwarp(Data{2,1},tform,'OutputView',imref2d(sizeim)); %move the green channel

% %% optional figure to check alignment
% figure
% imagesc(framesg_reg(:,:,1));axis image
% caxis([60 160])
% figure;imagesc(Data{2,2}(:,:,1)); axis image

% %% option to save transformation data if-needed
% filenamea = [fname,'_align.mat'];
% save([fpath,filenamea],'movingPoints','fixedPoints','tform');

%% 3. ID and track where peroxisomes are
% Using troika SPT code from Shuang, et al. PCCP 2014 
% GFP channel has peroxisomes - only track there
% LK20211215 - try to just localize in summed frame, not with tracking

toTrack=sum(framesg_reg,3);
toTrackBckCorrect=toTrack-mean(mean(toTrack));

% parameters for Troika stored in variable "e"
e.start_frame = 1; % the first frame to be analyzed
e.Gauss_width = 3;% estimate the Gaussian width of your PSF, not need to be integer
%the fitting region will be 4*Gauss_width+1 pixels
e.wide = 2;% the wide threshold used to define a local maximum
e.local_thd = true;% true if want to calculate local threshold, false if global threshold is fine.

% identify particles for each frame
im=toTrackBckCorrect;
position(1).p = particle_identify_highTresh(im, e.local_thd, e.Gauss_width, e.wide); %set threshold to 8x instead of 3x

trjR=[];
for p = 1 : numel(position(1).p(:,1))
    trjR(1, :, p) = position(1).p(p, :);
end
trjR_noMove=trjR;

% %% Optional code to uncomment if you would like to plot peroxisome
% %% locations over the raw data image to check ID-ing

% %% View found particles
% for i=1:150
%     imagesc(framesg_reg(:,:,i))
%     axis image
%     hold on
%     for j=1:size(trjR,3)
%         plot(trjR(1,1,j),trjR(1,2,j),'ro')
%     end
%     title(num2str(i))
%     caxis([100 150])
%     pause(0.01)
%     hold off
% end

% %% Map trajectories
% disp('mapping')
% % initial the trajectory matrix: first diminsion is time, second is x, y,
% % and Gaussian width, third dimension is the particle number.
% for p = 1 : numel(position(1).p(:,1))
%     trjR(1, :, p) = position(1).p(p, :);
% end
% map_new = [1: p; 1: p]';% match the map between two frames to particle numbers in trjR
% termied = [];% record the terminated particles
% for time = 1 : numel(position)-1
%     time;
%     map1_2 = mapping_frames(position(time).p, position(time+1).p, e.search_r);
%     [map_new, termied] = linking_map(map_new, map1_2, position, time, termied, e.search_r);
%     trjR = map2trj(trjR, map_new, position(time + 1).p, time + 1);
% end

%% 4. Run cross correlation on areas with peroxisomes
for i=1:size(trjR_noMove,3)
    peroxisomeLocatxy=[round(mean(trjR_noMove(:,1,i))),round(mean(trjR_noMove(:,2,i)))];  %find the average central pixel
    rangex=peroxisomeLocatxy(1)-8:peroxisomeLocatxy(1)+8; %we're interested in the 17 x 17 pixel area around the center
    rangey=peroxisomeLocatxy(2)-8:peroxisomeLocatxy(2)+8;

    if rangex(1)<1 || rangex(17)>size(Data{2,2},2)||rangey(1)<1 || rangey(17)>size(Data{2,2},2)
        XCresults(1,i)={[]};
        XCresults(2,i)={[]};
        XCresults(3,i)={[]};
    else
        for j=1:17
            for k=1:17    
                greendata=framesg_reg(rangey(j),rangex(k),:);
                greendata=reshape(greendata,1,size(greendata,3));
                greendata_corr=greendata-mean(greendata);
                reddata=Data{2,2}(rangey(j),rangex(k),:);
                reddata=reshape(reddata,1,size(reddata,3));
                reddata_corr=reddata-mean(reddata);
                XCdata_raw=xcorr(greendata_corr',reddata_corr','coeff'); %key step! This is where the correlation is performed
                XCdata=XCdata_raw(size(greendata,2):size(XCdata_raw,1));
                allXCdata{j,k}=XCdata;
                XClag1(j,k)=XCdata(1,1); %since XC, not AC, can use lag at tau=0
            end
        end
        % %% XCresults stores the important output data as a [3 x number of
        % %% peroxosiomes] cell
        XCresults(1,i)={peroxisomeLocatxy}; % 2x1 vector the center x,y position of the peroxisome in registered GFP data
        XCresults(2,i)={allXCdata}; %17x17 cell that indicates the pixel around the center position of the peroxisome; each entry is 150x1 G(tau) XC data
        XCresults(3,i)={XClag1}; %17x17 matrix of the magnitude of G(0) at each pixel in the peroxisome
    end

end

%% 5. Save data
if saveyn==1
    savename=[fname,'Analyzed.mat'];
    save(strcat(fpath,savename),'tform','XCresults','trjR_noMove','XCresults');
    % saves XCresults, along with the transformation for the image
    % registration of the two channels and the peroxisome locations from
    % the single particle tracking
end

%% 6. Optional scripts for making figures/post analysis of the outputted data
% %% Figure of all identified peroxisomes at G_X_C(0)
% figure
% nperox = size(XCresults,2);
% nrows=ceil(nperox/6);
% 
% for i=1:nperox
%     subplot(nrows,6,i)
%     imagesc(XCresults{3,i})%./max(max(XCresults{3,i})))
%     axis image; xticks(''); yticks(''); title([num2str(XCresults{1,i}),', ',num2str(i)]);
%     colormap(jet)
%     caxis([-0.1 0.6])
% end

% %% Movie with GFP and mCherry channels shown too
% figure
% indexPerox=3; %peroxisome you want to make a movie of
% for i=1:150
%     toPlot=zeros(17,17);
%     for j=1:17
%         for k=1:17
%             toPlot(j,k)=XCresults{2,indexPerox}{j,k}(i,1);
%         end
%     end    
%     imagesc(toPlot./max(max(XCresults{3,indexPerox})))
%     axis image; xticks(''); yticks('');
%     caxis([0 1])
%     tsec=i*0.15;
%     title(strcat(num2str(tsec),' s'))
%     pause(0.1)
%     hold off
% end

% %% Plot G_XC vs tau for a given pixel in a given peroxisome
% figure
% hold on
% indexPerox=1; %peroxisome of interest
% indexPixelrow=10; %pixel of interest
% indexPixelcol=10;
% times=[1:150]*0.15;
% plot(times,XCresults{2,indexPerox}{indexPixelrow,indexPixelcol},'LineWidth',2)
% set(gca,'LineWidth',1.5,'FontName','Malgun Gothic','FontSize',14)
% xlabel('\tau (s)')
% ylabel('G(\tau)')
% ylim([-0.4 0.6])
% xlim([0 22.5])

% %% Calculate the number peroxisomes with high cross correlation >0.5 at
% %% tau=0
% nperox=size(XCresults,2)
% nhigh=0;
% for i=1:size(XCresults,2)
%     if max(max(XCresults{3,i}))>0.5
%         
%         nhigh=nhigh+1;
%     end  
% end
% nhigh
% perchigh=nhigh./size(XCresults,2)
