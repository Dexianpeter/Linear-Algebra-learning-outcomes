%
% Homework 6 Matlab assignment
%
% Example code
% Image registration
% by linear transformations
%
% Linear Algebra
% NCKU, CSIE 
%
% Ming-Long Wu, Ph.D.
% 2024/6
%

clear all;
close all;
clc;


%%%%%%%
%
% 1.1 Loading images, key points, and preprocessing
%
%%%%%%%

% read images
IMA_1 = imread('./Data/A01_1.jpg');
IMA_2 = imread('./Data/A01_2.jpg');

% resize to 0.25x
IMA_1 = imresize(IMA_1, 0.25);
IMA_2 = imresize(IMA_2, 0.25);

% read labeled points
fid1 = fopen('./Data/control_points_A01_1_2.txt');
pts = fscanf(fid1,'%f');
fclose(fid1);

% shrink coordinates also by 0.25
pts = round(pts/4);

% sort X, Y coordinates of key points
X1 = pts(1:4:end);
Y1 = pts(2:4:end);
X2 = pts(3:4:end);
Y2 = pts(4:4:end);

% make a copy of images to add green key points
IMA_1m = IMA_1;
IMA_2m = IMA_2;

for i = 1 : length(X1)

    % each key point point is marked as a 11x11 green region for visualization
    % mark points for IMA_1m
    IMA_1m(Y1(i)-5:Y1(i)+5,X1(i)-5:X1(i)+5,:) = repmat(reshape([0, 255, 0], [1, 1, 3]),[11 11 1]);
    % mark points for IMA_2m
    IMA_2m(Y2(i)-5:Y2(i)+5,X2(i)-5:X2(i)+5,:) = repmat(reshape([0, 255, 0], [1, 1, 3]),[11 11 1]);
    
end

% display images to check if label points are correctly marked
figure(100);
% diplay A01-1.jpg
subplot(2,3,1);
imagesc(IMA_1m);
title('A01-1.jpg');
grid on;
xticks([0:120:728]);
yticks([0:120:728]);
set(gca,'dataaspectratio', [1 1 1],'linewidth', 2, 'GridColor', [1 1 1], 'GridAlpha', 0.7, 'FontSize', 18);
% display A01-2.jpg
subplot(2,3,2);
imagesc(IMA_2m);
title('A01-2.jpg');
grid on;
xticks([0:120:728]);
yticks([0:120:728]);
set(gca,'dataaspectratio', [1 1 1],'linewidth', 2, 'GridColor', [1 1 1], 'GridAlpha', 0.7, 'FontSize', 18);


%%%%%%%
%
% 1.2 Finding affine transformation coefficients
%
%%%%%%%

% Solve two least square problems here.

N = length(X1);
A = [X2, Y2, ones(N,1)];
b1 = X1;
T_tmp1 = lsqr(A, b1, 1e-6, 100);
b2 = Y1;
T_tmp2 = lsqr(A, b2, 1e-6, 100);

T = [T_tmp1'; T_tmp2'; 0, 0, 1];
fprintf("Affine transformation = \n");
disp(T);

% Calculate root mean square error
XY2_before = [X2'; Y2';ones(1,N)];
XY2_after = T * XY2_before;

err = 0;
for i = 1 : N
    err = sqrt( (XY2_after(1, i)- X1(i))^2 + (XY2_after(2, i)- Y1(i))^2 ) + err;
end
rmse = sqrt(err/N);
fprintf("rmse = %f (pixels)\n", rmse);
fprintf("(2)\n");
fprintf("把圖形旋轉後，誤差已經小於1個像素了，\n");
fprintf("可以表現出準確度非常的高。\n");
%%%%%%%
%
% 1.3 Applying affine transformation to A01-2.jpg
%
%%%%%%%

[H, W, ~] = size(IMA_2m);
[X2_generate, Y2_generate] = meshgrid(0:W-1, 0:H-1);
X2_before = X2_generate(:);
Y2_before = Y2_generate(:);

% Apply affine transformation.
XY2_before = [X2_before'; Y2_before'; ones(1, numel(X2_before))];
XY2_after = T *XY2_before;
[X_new, Y_new] = meshgrid(0:W-1, 0:H-1); 

I2_double = im2double(IMA_2m);
R = I2_double(:,:,1);  
r2 = R(:); 
G = I2_double(:,:,2); 
g2 = G(:);
B = I2_double(:,:,3); 
b2 = B(:);

R_new = griddata(XY2_after(1,:), XY2_after(2,:), r2, X_new, Y_new, 'linear');
G_new = griddata(XY2_after(1,:), XY2_after(2,:), g2, X_new, Y_new, 'linear');
B_new = griddata(XY2_after(1,:), XY2_after(2,:), b2, X_new, Y_new, 'linear');

IMA2_reg = cat(3, R_new, G_new, B_new);

%%%%%%%
%
% 1.4 Check image registration
%
%%%%%%%

% Display whole images and cropped images to check image registration.
subplot(2,3,3);
imagesc(IMA2_reg);
title('A01-2-registration.jpg');
grid on;
xticks([0:120:728]);
yticks([0:120:728]);
set(gca,'dataaspectratio', [1 1 1],'linewidth', 2, 'GridColor', [1 1 1], 'GridAlpha', 0.7, 'FontSize', 18);

ROI_x = 300:300+363;
ROI_y = 200:200+363;
I1_ROI = IMA_1m(ROI_y, ROI_x, :);  
I2_ROI = IMA_2m(ROI_y, ROI_x, :);    
I2reg_ROI = IMA2_reg(ROI_y, ROI_x, :);

subplot(2,3,4);
imagesc(I1_ROI);
grid on;
title('A01-1–cropped');
set(gca,'dataaspectratio', [1 1 1],'linewidth', 2, 'GridColor', [1 1 1], 'GridAlpha', 0.7, 'FontSize', 18);

subplot(2,3,5);
imagesc(I2_ROI);
grid on;
title('A01-2–cropped');
set(gca,'dataaspectratio', [1 1 1],'linewidth', 2, 'GridColor', [1 1 1], 'GridAlpha', 0.7, 'FontSize', 18);

subplot(2,3,6);
imagesc(I2reg_ROI);
grid on;
title('A01-2–registration-cropped');
set(gca,'dataaspectratio', [1 1 1],'linewidth', 2, 'GridColor', [1 1 1], 'GridAlpha', 0.7, 'FontSize', 18);

fprintf("(3)\n")
fprintf("從A01-2cropped和A01-1中可以看到，\n" + ...
    "旋轉後綠色的點和A01中綠色點相當接近，\n" + ...
    "以血管的大小來說，這樣的誤差範圍應是能接受的。\n" +...
    "所以我認為這是非常有效的做法。\n");