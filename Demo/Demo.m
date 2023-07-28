%% Demo using stock image
% Images
ref = uint16(imread("cameraman.tif"));

% Stack
nframes = 10;
stack = zeros(size(ref,1), size(ref,2), nframes);

% Shearing and cropping
for i = 1 : nframes
    % Shearing
    sh = randn*0.1;
    shear = [1 0 0; sh 1 0; 0 0 1]; 
    tform = maketform('affine',shear);
    moving = imtransform(ref,tform); 

    % Cropping
    crop = round(size(moving,2) - size(ref,2));
    cropv = [round(crop/2), crop-round(crop/2)];
    moving = moving(:, 1+cropv(1):end-cropv(2));
    stack(:,:,i) = moving;
end

% Show
figure
subplot(1,2,1)
imshow(ref,[]);
title('Reference')
subplot(1,2,2)
imshow(mean(stack,3),[]);
title('Zmean')

%% Register
[stack_reg, ~] = DemonsRegOneshot(stack, 'ref', ref, 'parfor', false, 'toseg', false);

% Show
figure('Position', [50 100 1200 500]);
subplot(1,3,1)
imshow(ref,[]);
title('Reference')

subplot(1,3,2)
imshow(mean(stack,3),[]);
title('Pre-registration zmean')

subplot(1,3,3)
imshow(mean(stack_reg,3),[]);
title('Post-registration zmean')