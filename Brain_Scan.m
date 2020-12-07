clear       %no variables
close all   %no figures
clc         %empty command window
%%
%the filename convention 

prompt = 'Type Tumour Folder as (Tumour-#/1-):';
prefix = input(prompt, 's');
ext = '.dcm';
fnum = input('Enter the Number of Brain Slices (write vector in square brackets:');

%first files name in the dataset 
fname = [prefix num2str(fnum(1)) ext];

%obtain and examine the fileheader 
info = dicominfo(fname);

%extract the size information 
voxel_size = [info.PixelSpacing; info.SliceThickness];

%read the image and populate an XYZ vector space 
hWaitBar = waitbar(0,'Reading DICOM files');
for i=length(fnum):-1:1
  fname = [prefix num2str(fnum(i)) ext];
  D(:,:,i) = uint16(dicomread(fname));
  waitbar((length(fnum)-i)/length(fnum))
end
delete(hWaitBar);

%imshow(fname)

%%
%data visualization
%explore image data using Image Viewer
i = input('Which MRI slice do you want to analyze:');  
im = squeeze(D(:,:,i));
max_level = double(max(D(:))); 
figure(1)
imshow(im, [0 max_level]);

%%

%convert to grayscale
gray_im = mat2gray(im);
%adjusted grayscale

adj_im = imadjust(gray_im);
figure(2)
imshow(adj_im);
title('Enhanced Grayscale MRI');

%%

T1 = graythresh(adj_im);
%Highpass Filter 
% Getting Fourier Transform of the input_image 
% using MATLAB library function fft2 (2D fast fourier transform)   
[M, N] = size(adj_im);
FT_img = fft2(double(adj_im)); 
  
% Assign Cut-off Frequency   
D0 = T1+0.77; % one can change this value accordingly 
  
% Designing filter 
u = 0:(M-1); 
idx = find(u>M/2); 
u(idx) = u(idx)-M; 
v = 0:(N-1); 
idy = find(v>N/2); 
v(idy) = v(idy)-N; 
  
% Matrix V with each row is a copy of v, and matrix U  
% with each column is a copy of u 
[V, U] = meshgrid(v, u); 
  
% Calculating Euclidean Distance 
D = sqrt(U.^2+V.^2); 
  
% Comparing with the cut-off frequency and  
% determining the filtering mask 
H = double(D > D0); 
  
% Convolution between the Fourier Transformed image and the mask 
G = H.*FT_img; 
  
% Getting the resultant image by Inverse Fourier Transform 
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)   
output_im = real(ifft2(double(G))); 

% Displaying Input Image and Output Image  
figure(3)
imshow(output_im); 
%%
%MEDIAN Filter
%convert adj_m to double
Med_MRI = medfilt2(output_im, [2 2]);
figure(4)
hold on
imshow(Med_MRI)
hold off

%determine minimum threshold 
T2 = graythresh(Med_MRI);

%Binarize to obtain tru B/W
im_binarize = imbinarize(Med_MRI, T1);
figure(5)
imshow(im_binarize);
title('Segmented MRI');
%%
%morphology 
im_morph = bwmorph(im_binarize, 'bridge', inf);
imshow(im_morph);
%%
