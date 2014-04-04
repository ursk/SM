

% make a really cool colormap:
r=[ones(1,128),        linspace(1,0, 128)];
g=[linspace(0,1, 128), linspace(1,0, 128)];
b=[linspace(0,1, 128), ones(1,128)       ];

map=[r' g' b'];
colormap(map)
%map(1,:)=0; %black extrema
%map(256,:)=0;