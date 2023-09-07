mov = VideoReader('SundialTimeLapse720Short.mp4');   %# use mmreader on older versions
for i=1:mov.NumberOfFrames
    img = read(mov,i);
    imwrite(img, strcat('Sundial720/image',num2str(i),'.jpg'));
end

% videoFReader = vision.VideoFileReader('chairTimelapse720.mp4','VideoOutputDataType', 'uint8');
% videoPlayer = vision.VideoPlayer;
% frame = step(videoFReader);
% a = energyRGB(frame);
% for i=1:100
%    frame(i,i,:) = [0,0,255]; 
% end
% 
% 
% imshow(a);


% while ~isDone(videoFReader)
%    frame = step(videoFReader);
%    ind2rgb(X,map) 
% %    disp(frame(100,100,:));
% %    energyWindow = energyRGB(frame);
% %    step(videoPlayer,frame);
% end

% release(videoFReader);
% release(videoPlayer);

%% This is the preprocessing stage for calculating the shadow line angles in the first image

%This part will be automized later by getting information from ground
%shadow detection code

% firstImage = imread('Sundial/image1','jpg');
% trackingWindow = firstImage(430:559,575:874,1:3);
% energyWindow = energyRGB(trackingWindow);
% %imshow(energyWindow);
% 
% firstStart = [1,126];
% firstEnd = [130,261];
% 
% secondStart = [1,247];
% secondEnd = [44,300];
% 
% firstAngle = atand((firstEnd(1,1)-firstStart(1,1)+1)/(firstEnd(1,2)-firstStart(1,2)+1));
% secondAngle = atand((secondEnd(1,1)-secondStart(1,1)+1)/(secondEnd(1,2)-secondStart(1,2)+1));
% initialAverageAngle = (firstAngle + secondAngle)/2;
% 
% %% This part is to read and evaluate the upcoming images
% 
% tic
% difference = zeros(58,1);
% for i=2:58
%     
%     nextImage = imread(strcat('Sundial/image', num2str(i),'.jpg'));
%     trackingWindow = nextImage(430:559,575:874,1:3);
%     energyWindow = energyRGB(trackingWindow);
%     %imshow(energyWindow);
% 
%     firstStart = [1, find(energyWindow(1,:) == max(energyWindow(1, (firstStart(1,2)-10):(firstStart(1,2)+10))),1)];
%     firstEnd = [130, find(energyWindow(130,:) == max(energyWindow(130, (firstEnd(1,2)-10):(firstEnd(1,2)+10))),1)];
% 
%     secondStart = [1, find(energyWindow(1,:) == max(energyWindow(1, (secondStart(1,2)-10):(secondStart(1,2)+10))),1)];
%     secondEnd = [find(energyWindow(:,300) == max(energyWindow((secondEnd(1,1)-10):(secondEnd(1,1)+10), 300)),1), 300];
% 
%     firstAngle = atand((firstEnd(1,1)-firstStart(1,1)+1)/(firstEnd(1,2)-firstStart(1,2)+1));
%     secondAngle = atand((secondEnd(1,1)-secondStart(1,1)+1)/(secondEnd(1,2)-secondStart(1,2)+1));
%     currentAverageAngle = (firstAngle + secondAngle)/2;
% 
%     difference(i) = currentAverageAngle - initialAverageAngle;
%     initialAverageAngle = currentAverageAngle;
% 
% end
% b = sum(difference);
% toc