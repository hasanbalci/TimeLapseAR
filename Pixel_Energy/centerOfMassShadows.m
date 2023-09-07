tic
firstImage = imread('Chair/image1','jpg');
energyImage1 = energyRGB(firstImage(560:719,700:1034,:));
% cenOfMass = centerOfMass(energyImage);
% imshow(energyImage1);

% imshow(energyImage2);

cenOfMass = zeros(30,2);
for i=2:2:600
    
    nextImage = imread(strcat('Chair/image', num2str(i),'.jpg'));
    energyImage2 = energyRGB(nextImage(560:719,700:1034,:));
    diffImage = energyImage2 - energyImage1;
    diffImage(diffImage(:,:)<10) = 0;
    cenOfMass(floor(i/2)+1,:)  = centerOfMass(diffImage);
    energyImage1 = energyImage2;
    
end
% imshow(diffImage);

figure; imagesc(diffImage); colormap gray; axis image
hold on; plot(cenOfMass(:,2),cenOfMass(:,1),'rx');
toc