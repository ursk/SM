% General data for all models
% 
% Arguments
% - samples: number of image patches to sample, e.g. 10000
% - winsize: sampling window, e.g. 8 pixels
% - rdim: PCA dimensionality reduction. If rdim=winsize^2, ZPW is used
% - cgc: neighborhood size for gain control. If 0, no gain control is performed
%
% Return va


function [X, wM, dwM] = data( samples, winsize, rdim, cgc, seed);

t=round(samples*1);
if nargin==5; rand('seed', seed); else; rand('seed', 0);  end;
N=256; %square image size

%----------------------------------------------------------------------
% 1. Sample patches from images
%----------------------------------------------------------------------

dataNum = 13;% We have a total of 13 images.
getsample = floor(t/dataNum);% This is how many patches to take per image
X = zeros(winsize^2,getsample*dataNum);% Initialize the matrix to hold the patches


sampleNum = 1;  
for i=(1:dataNum) % loop for images
    if i==dataNum, getsample = t-sampleNum+1; end  % Even things out (take enough from last image)
  % Load the image
    I = imread(['./data/' num2str(i) '.tiff']);  
    I=I(1:N,1:N);
    I = double(I);
    I = I-mean(mean(I));
    I = I/sqrt(mean(mean(I.^2)));
    
  if cgc
      % whiten the whole image for CGC
      imf=fftshift(fft2(I)); 
      [fx fy] = meshgrid(-N/2:N/2-1);% generate the filter with highpass roll-off
      [theta rho]=cart2pol(fx,fy);
      filtf=rho.*exp(-.5*(rho/(.7*N/2)).^2); 
      imwf=filtf.*imf;  %apply filter in frequency domain
      imw=real(ifft2(fftshift(imwf)));% go back to space domain
      imw = imw-mean(mean(imw)); % zero mean
      imw = imw/sqrt(mean(mean(imw.^2))); % unit variance
   
      %perform gain control
      [x y]= meshgrid(-cgc/2:cgc/2-1);
      G=exp(-.5*((x.^2+y.^2)/(cgc/2)^2));
      G=G/sum(G(:));
      imv=conv2(imw.^2,G,'same'); % convolve image with kernel to find contrast
      I=imw./sqrt(imv+.1);
  end
  
  sizex = size(I,2); sizey = size(I,1);
   
  %take samples from the image
  sizex = size(I,2); sizey = size(I,1);
  posx = floor(rand(1,getsample)*(sizex-winsize-2))+1;
  posy = floor(rand(1,getsample)*(sizey-winsize-1))+1;
  
  for j=1:getsample %patches per image
    X(:,sampleNum) = reshape( I(posy(1,j):posy(1,j)+winsize-1, ...
			posx(1,j):posx(1,j)+winsize-1),[winsize^2 1]);
    sampleNum=sampleNum+1;
  end 
  
end % going though images



%----------------------------------------------------------------------
% 2. Whiten
%----------------------------------------------------------------------
if rdim < winsize^2; X = X - repmat(mean(X), winsize^2, 1); end % remove DC for PCA whitening

[E, D] = eig(cov(X')); % Eigenvalue decomposition

[dummy,order] = sort(diag(-D));
% E and D are calculated once for rDim, then used twice
E = E(:,order(1:rdim));
d = diag(D); 
d = real(d.^(-0.5));
D = diag(d(order(1:rdim)));

if rdim == winsize^2
  wM =E*D*E'; % zero phase whitening
  dwM=inv(wM);
else 
  wM = D*E';
  dwM = E * inv(D);
end

X= wM*X; %whiten the data

return;

