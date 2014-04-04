function visual( A, mag, cols, norma )
% visual - display a basis for image patches
% Based on Patrik Hoyer's ImageICA package
%
% A        the basis, with patches as column vectors
% mag      magnification factor
% cols     number of columns (x-dimension of map)
% norma    normalize each patch to cove the full scale

% Get maximum absolute value (it represents white or black; zero is gray)
maxi=max(max((A)));
mini=min(min((A)));

% This is the side of the window
dim = sqrt(size(A,1));

% Helpful quantities
dimm = dim-1;
dimp = dim+1;
rows = size(A,2)/cols;
if rows-floor(rows)~=0, error('Fractional number of rows!'); end

% Initialization of the image - WAS MAXI * ...
I = maxi*ones(dim*rows+rows-1,dim*cols+cols-1); 

for i=0:rows-1
  for j=0:cols-1
    
    % This sets the patch
    patch=reshape(A(:,i*cols+j+1),[dim dim]);
    if nargin==4; patch=patch/ max(abs(patch(:))); end % normalize each patch
    I(i*dimp+1:i*dimp+dim,j*dimp+1:j*dimp+dim) = patch	;
  end
end

Ir = imresize(I+1,mag,'nearest')-1;

%figure;
%colormap(gray(256));
%iptsetpref('ImshowBorder','tight'); 
%subplot('position',[0,0,1,1]);
%maxi=max((I(:)));
%mini=min((I(:)));
imagesc(Ir); axis image 
%colormap(jet);
french
%colorbar
truesize;  
%colorbar
