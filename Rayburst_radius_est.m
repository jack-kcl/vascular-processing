%% main function: reads in the skeleton (and Analyze skeleton data), detects radii

% folder leads to the mask, with skel underneath it
function [R pt segs] = Rayburst_radius_est(folder)

%numrays = 60;
subdiv = 2;
trunc = 0.5; %  the delete-axis should be aligned with vessel axis!

if 0 % cylinder phantom
    mask = zeros(101,101,101);
    [X Y Z] = meshgrid(-50:50);
    R= folder;  % using "folder" as cylinder radius
    mask((Y-Z).^2+(Z-X).^2+(X-Y).^2<=3*R^2) = 255;
    pt=[51 51 49];
end
if 1 % use real data
    % -- import the image stack
    fprintf('importing skeleton...\n');
    skel = import_stack([folder '\skel']);
    fprintf('importing vessel mask...\n');
    mask = import_stack(folder);
    mask = double(mask);
end

if 0
    % -- define points for measurement
    fprintf('detecting centreline points...\n');
    allnodes = find(skel);
    pt = zeros(size(allnodes,1),3);
    [pt(:,1),pt(:,2),pt(:,3)] = ind2sub(size(skel),allnodes);
end
if 1
    % --- use skeleton import properly
    fprintf('importing skeleton centreline points...\n');  % from Fiji
    if 0  % import Fiji AnalyseSkeleton 2D/3D .xlsx file
        [knots,nodes,segs] = import_skeleton([folder filesep 'skel']);
        segs=knots;  
    end
    if 1  % reconstruct segs directly from skel image stack
        segs = AnalyseSkeleton([folder filesep 'skel'], 'skel');
    end
    ns=length(segs); npt=0; 
    for j=1:ns
        npt=npt+size(segs{j},1);
    end
    % make pt array (can reduce, if necessary)
    R = zeros(npt,1); pt=zeros(npt,3); jj=0;
    for j=1:ns
        pt(jj+1:jj+size(segs{j},1),:) = segs{j};
        jj=jj+size(segs{j},1);
    end
end

% precalc sampling rays
%rays = sampling_rays(numrays);  % alternative samplers
rays = sampling_rays2(subdiv);
quart=round(length(rays)/4);
%-- probably want to adjust by R2 = (R-0.42)/1.1 (in voxels)

% estimate radius: for each direction, return two rad measures (+/-)
fprintf('measuring radii...\n'); t0=tic;
for c = 1:size(pt,1)
    fprintf('processing point %i of %i...\n',c,size(pt,1));
    [rads rf rb] = measure(mask, pt(c,:), rays);
    % how to pick the true radius?
    r=sort(rads);   R(c)=r(quart);
end
t1=toc(t0); fprintf('%g seconds elapsed\n',t1);

% export for viewing
vars{1}.name='radius';
vars{1}.data=R;
vars{1}.ncomp=1;
here=pwd; cd(folder);
write_exnode_file('skelrads',pt,length(pt),3,vars,1,[1]);
write_vtk_file('skelrads',length(pt),0,'tri',1,pt,[], [1],vars);
cd(here);
end
%% --- measure radii
function [rad rf rb] = measure(im, pt, rays)
    % it's assumed that rays only contains half the points
    nr = size(rays,1);
    rf = zeros(nr,1); rb = rf; rad = rf;
    % terminating condition: first to hit 0.5
    % some will not hit it (axially aligned)
    % adaptive detection needed
    thresh = 0.5;
    for i=1:nr
        % do 10 at a time?
        scales = [1:.25:12]';
        ray(:,1) = scales*rays(i,1);
        ray(:,2) = scales*rays(i,2);
        ray(:,3) = scales*rays(i,3);
        intf = mirt3D_mexinterp(im, ray(:,2)+pt(2), ray(:,1)+pt(1), ray(:,3)+pt(3));
        intb = mirt3D_mexinterp(im, -ray(:,2)+pt(2), -ray(:,1)+pt(1), -ray(:,3)+pt(3));
        % detect crossing
        try
        rf(i) = min(find(intf<thresh));
        rb(i) = min(find(intb<thresh));
        rad(i) = (scales(rf(i))+scales(rb(i)))/2;
        catch
            rad(i) = max(scales);
        end
    end
end
%% --- golden spiral 
function p = sampling_rays(numrays)
p=zeros(3,numrays);
inc = pi/(3.0-sqrt(5.0));
off = 1.0/(2.0*numrays) - 1;
for k=0:numrays-1
    z = k/numrays + off;
    r = sqrt(1.0-z*z);
    phi = k*inc;
    p(:,k+1) = [cos(phi)*r;
              sin(phi)*r;
              z];
end
% only fills half the sphere.. perfect
end
%% --- dodecahedron
function p = sampling_rays2(subdiv)
phi = (1+sqrt(5))/2;
b = 1 / phi ; 
c = 2 - phi ;
V=[-1 -c 0
   -1 c 0
   -b -b -b
   -b -b b
   -b b -b
   -b b b
   -c 0 -1
   -c 0 1
   0 -1 -c
   0 -1 c
   0 1 -c
   0 1 c
   c 0 -1
   c 0 1
   b -b -b
   b -b b
   b b -b
   b b b
   1 -c 0
   1 c 0];
F=[14    18    12     6     8
     8     4    10    16    14
    13    15     9     3     7
     7     5    11    17    13
    11    12    18    20    17
    12    11     5     2     6
     9    10     4     1     3
    10     9    15    19    16
    20    18    14    16    19
    19    15    13    17    20
     2     5     7     3     1
     1     4     8     6     2];
% have to make triangles
nv=size(V,1); nf=size(F,1);
V(nv+nf,:)=[0 0 0];
for i=1:nf
    V(nv+i,:) = mean(V(F(i,:),:));
end
p = normalise(V);
% make elements
el=[]; ne=0;
for i=1:nf
    for j=1:4
        ne=ne+1;
        el(ne,:)=[F(i,j) F(i,j+1) i+nv];
    end
    ne=ne+1;
    el(ne,:) = [F(i,5) F(i,1) i+nv];
end
% subdivide now
for i=1:subdiv  % 1:122 points, 2:482 points, 3:1922 points.
    [p el] = subdivide_tri(p, el);
    p = normalise(p);
end
% remove half of them
p(p(:,3)<0,:)=[];
p(p(:,3)<1e-5 & p(:,2)<0,:)=[];
p(p(:,3)<1e-5 & p(:,2)<1e-5 & p(:,1)<0,:)=[];
%plot3(rays(:,1),rays(:,2),rays(:,3),'.'); axis equal
%hold on; plot3(-rays(:,1),-rays(:,2),-rays(:,3),'rx');
end
%% --- normalise all coordinates
function V = normalise(V)
r = sqrt(V(:,1).^2 + V(:,2).^2 + V(:,3).^2);
V(:,1) = V(:,1)./r;
V(:,2) = V(:,2)./r;
V(:,3) = V(:,3)./r;
end
%% --- 
function im = import_stack(folder)
here = pwd;
cd(folder);
try
    ls = dir('*.png');
    n = length(ls);
    info = imfinfo(ls(1).name);
    switch info(1).BitDepth
    case {8}
        type='uint8';
    case {16}
        type='uint16';
    case {32}
        type='double';
    otherwise
        error('cannot determine image bit depth');
    end
    im = zeros(info(1).Height, info(1).Width, n, type);
    for i=1:n
        im(:,:,i) = imread(ls(i).name, info(1).Format);
    end
catch
    disp('could not import images...');
end
cd(here);
end

%% ---
function [v1 f1] = subdivide_tri( xyz, faces )
%
% Vectorized Triangle Subdivision: split each triangle 
% input face into four new triangles. 
% 
% usage: [v1 f1] = subdivide( xyz, faces );
%
%  author:  Peter A. Karasev     25 Nov 2009


numverts = size(xyz,1);
numfaces = size(faces,1);
disp(['Input mesh: ' num2str(numfaces) ' triangles, ' ... 
    num2str(numverts) ' vertices.']);

fk1 = faces(:,1);
fk2 = faces(:,2);
fk3 = faces(:,3);

% create averages of pairs of vertices (k1,k2), (k2,k3), (k3,k1)
    m1x = (xyz( fk1,1) + xyz( fk2,1) )/2;
    m1y = (xyz( fk1,2) + xyz( fk2,2) )/2;
    m1z = (xyz( fk1,3) + xyz( fk2,3) )/2;
    
    m2x = (xyz( fk2,1) + xyz( fk3,1) )/2;
    m2y = (xyz( fk2,2) + xyz( fk3,2) )/2;
    m2z = (xyz( fk2,3) + xyz( fk3,3) )/2;
    
    m3x = (xyz( fk3,1) + xyz( fk1,1) )/2;
    m3y = (xyz( fk3,2) + xyz( fk1,2) )/2;
    m3z = (xyz( fk3,3) + xyz( fk1,3) )/2;

    
vnew = [ [m1x m1y m1z]; [m2x m2y m2z]; [m3x m3y m3z] ];
clear m1x m1y m1z m2x m2y m2z m3x m3y m3z
[vnew_ ii jj] = unique(vnew, 'rows' );

clear vnew; 
m1 = jj(1:numfaces)+numverts;
m2 = jj(numfaces+1:2*numfaces)+numverts;
m3 = jj(2*numfaces+1:3*numfaces)+numverts;

tri1 = [fk1 m1 m3];
tri2 = [fk2 m2 m1];
tri3 = [ m1 m2 m3];
tri4 = [m2 fk3 m3];
clear m1 m2 m3 fk1 fk2 fk3
 
v1 = [xyz; vnew_]; % the new vertices
f1 = [tri1; tri2; tri3; tri4]; % the new faces
disp(['Output mesh: ' num2str(size(f1,1)) ' triangles, ' ... 
    num2str(size(v1,1))  ' vertices.']);
end
