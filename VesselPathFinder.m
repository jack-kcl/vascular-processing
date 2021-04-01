function [knots,nodes,segs] = VesselPathFinder(sknodes,sksegs,skelvol)

skelvolp = zeros(size(skelvol,1)+2,size(skelvol,2)+2,size(skelvol,3)+2,'uint8');
skelvolp(2:size(skelvolp,1)-1,2:size(skelvolp,2)-1,2:size(skelvolp,3)-1) = skelvol;
clear skelvol;
skelvolp = logical(skelvolp); %skelvolp(skelvolp>0)=255;

% % ???
% skelvolp = logical(skelvol);

allnodes = find(skelvolp);
allnodesS = zeros(size(allnodes,1),3);
[allnodesS(:,1),allnodesS(:,2),allnodesS(:,3)] = ind2sub(size(skelvolp),allnodes);

%tic

xs = size(skelvolp,1); ys = size(skelvolp,2); zs = size(skelvolp,3);

allsegs = [];

for nodeind = 1:size(allnodesS,1)
    
    %disp(nodeind)
    
    clear n xb xf yb yf zb zf vol nze nzeS translate neighbors neighborsI ne nni
    
    n = allnodesS(nodeind,:); % node
    xb = n(1) - 1; %if xb < 1, xb = 1; end
    xf = n(1) + 1; %if xf > xs, xf = xs; end
    yb = n(2) - 1; %if yb < 1, yb = 1; end
    yf = n(2) + 1; %if yf > ys, yf = ys; end
    zb = n(3) - 1; %if zb < 1, zb = 1; end
    zf = n(3) + 1; %if zf > zs, zf = zs; end
    
    vol = skelvolp(xb:xf,yb:yf,zb:zf); % extract the 3x3x3 vol
    nze = find(vol>0); % find all non-zero elements in the vol
    nze(nze==14) = []; % remove the central voxel which is 14 in index <=> (2,2,2) in subscripts
if isempty(nze); continue; end
    [nzeS(:,1),nzeS(:,2),nzeS(:,3)] = ind2sub(size(vol),nze); % get subscripts of all non-zero elements
    translate = nzeS - repmat([2,2,2],size(nzeS,1),1); % repeat in rows dimension as many times as there are rows in nzeS, repeat in cols dimension only once
    neighbors = repmat(n,size(translate,1),1) + translate; % perform translation to get neighbor nodes
    neighborsI = sub2ind(size(skelvolp),neighbors(:,1),neighbors(:,2),neighbors(:,3));
    
    for ne = 1:size(neighborsI,1)
        nni = find(allnodes == neighborsI(ne));
        allsegs = [allsegs;nodeind,nni];
    end
    
end

disp('Finished creating the map network of the binary vasculature.')
%toc

%% Dijkstra Shortest Path

%tic

sknodesc = sknodes; clear sknodes
sknodes(:,1) = sknodesc(:,2); sknodes(:,2) = sknodesc(:,1); sknodes(:,3) = sknodesc(:,3); % flip x and y for imageJ -> matlab

sknodes(:,3) = sknodes(:,3) + 1; % add 1 in z dimension because 1 extra slice padded in beginning of z-dimension
sknodes(:,2) = sknodes(:,2) + 1; % add 1 in y dimension because 1 extra slice padded in beginning of y-dimension
sknodes(:,1) = sknodes(:,1) + 1; % add 1 in x dimension because 1 extra slice padded in beginning of x-dimension

sknodes(:,3) = sknodes(:,3) + 1; % add 1 in z dimension because imageJ starts counting at 0 and MATLAB starts counting at 1
sknodes(:,2) = sknodes(:,2) + 1; % add 1 in y dimension because imageJ starts counting at 0 and MATLAB starts counting at 1
sknodes(:,1) = sknodes(:,1) + 1; % add 1 in x dimension because imageJ starts counting at 0 and MATLAB starts counting at 1

clear sknodesc;

svs = size(skelvolp);
allnodesSID(:,1) = 1:size(allnodesS,1); allnodesSID(:,2:4) = allnodesS;
allsegsID(:,1) = 1:size(allsegs,1); allsegsID(:,2:3) = allsegs;

am = logical(sparse(size(allnodesSID,1),size(allnodesSID,1))); % sparse adjacency matrix
for seg = 1:size(allsegsID,1) % iterate over allsegs
    am(allsegsID(seg,2),allsegsID(seg,3)) = 1;
    am(allsegsID(seg,3),allsegsID(seg,2)) = 1;
    %disp(seg)
end
am = double(am);

sksegs_loops = logical(zeros(size(sksegs,1),1,'uint8'));
flag = 1;
for seg = 1:size(sksegs,1) % iterate over skeleton segs or rows in sksegs
    s = sksegs(seg,1);
    e = sksegs(seg,2);
    sn = sknodes(s,1:3);
    en = sknodes(e,1:3);
    %     sni = sub2ind(svs,sn(1),sn(2),sn(3));
    %     eni = sub2ind(svs,en(1),en(2),en(3));
    %     sniANI = find(allnodes == sni);
    %     eniANI = find(allnodes == eni);
    if all(sn == en)%((s == e)||all(sn == en)||(sni == eni)||(sniANI == eniANI))
        sksegs_loops(seg,1) = 1;
        if(flag == 1)
            disp('Found a looping seg');
            flag = 0;
        end
    end
    %disp(seg)
end
sksegs_orig = sksegs;
sksegs(sksegs_loops==1,:) = []; % remove every seg that is a node looping back onto itself

knots = cell(size(sksegs,1),1);
for seg = 1:size(sksegs,1) % i s
    clear s e sn en sni eni sniANI eniANI dist pathd pathg pathdm pred coords coordsIND
    s = sksegs(seg,1);
    e = sksegs(seg,2);
    sn = sknodes(s,1:3);
    en = sknodes(e,1:3);
    %sni = sub2ind(svs,sn(1),sn(2),sn(3));
    %eni = sub2ind(svs,en(1),en(2),en(3));
    %sniANI = find(allnodes == sni);
    %eniANI = find(allnodes == eni);
    sniANI = find((allnodesS(:,1)==sn(1))&(allnodesS(:,2)==sn(2))&(allnodesS(:,3)==sn(3)));
    eniANI = find((allnodesS(:,1)==en(1))&(allnodesS(:,2)==en(2))&(allnodesS(:,3)==en(3)));
    pathdm = shortestpathfinderBGL(am,sniANI,eniANI);
    coords = allnodesSID(pathdm,2:4); % column 1 contains nodeID, and columns 2:4 contain x,y,z coordinates
    coords = coords - 1; % subtract 1 from all three columns (dimensions) because 1 extra slice was added in each dimension for padding needed for skeleton paths tracing
    knots{seg} = coords;
    %disp(seg)
end

disp('Finished tracing each vessel segment.')

nodes = sknodes - 1; % subtract 1 from all three columns (dimensions) because 1 extra slice was added in each dimension for padding needed for skeleton paths tracing
segs = sksegs;
%toc

end