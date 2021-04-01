%% AnalyseSkeleton
%  The Fiji version removes loops, even though it's in the skel image!
%  Let's roll our own

function [segs juncs] = AnalyseSkeleton(folder, file)
t0=tic;

% -- read in the image stack
fprintf('importing the skeleton images...\n');
stack = import_png_stack(folder, file);
%stack(stack>0)=1;  % replace with 1
ind = find(stack);
[X Y Z] = ind2sub(size(stack),ind);

% -- find degrees of all chain nodes
fprintf('analysing nodes...\n');
[deg neigh] = find_deg_neigh(stack, X,Y,Z);

% -- mark the stack with the degrees
stack = mark_stack(stack, X,Y,Z, deg);

% -- find each segment
fprintf('tracing segments...\n');
[segs stack] = find_segs(stack, X,Y,Z, deg,neigh);

% -- also grab the junction nodes (why not)
ind = find(deg~=2);
juncs = [X(ind) Y(ind) Z(ind)];

t1=toc(t0);
fprintf('elapsed time: %g\n',t1);



% vars{1}.name='radius';
% vars{1}.data=zeros(length(X),1);
% vars{1}.ncomp=1;
% write_exnode_file('skelrads',[Y X Z],length(X),3,vars,1,[1]);

end

%% ----
function [deg neigh] = find_deg_neigh(stack, x,y,z)
n = length(x);
deg=zeros(n,1); neigh{n}=[];
for i=1:n
    vol = stack(x(i)-1:x(i)+1, y(i)-1:y(i)+1, z(i)-1:z(i)+1);  % assume we're not close to the image border
%     x0=max(x(i)-1, 1); y0=max(y(i)-1, 1); z0=max(z(i)-1, 1);
%     x1=min(x(i)+1, x(end)); y1=min(y(i)+1, y(end)); z1=min(z(i)+1, z(end));
%     vol = stack(x0:x1, y0:y1, z0:z1);
    vol(2,2,2) = 0;
    neigh{i} = find(vol);
    deg(i) = length(neigh{i});
end
end
%% ----
function stack = mark_stack(stack, X,Y,Z, deg)
for i=1:length(deg)
    stack(X(i),Y(i),Z(i))=deg(i);
end
end
%% ----
function [segs stack] = find_segs(stack, x,y,z, deg,neigh)
n = length(deg);
segs{1}=[]; ns=0;
% mapping between sub and ind
ind=[1:27];
[I J K]=ind2sub([3 3 3],ind); I=I-2; J=J-2; K=K-2; 
dist=abs(I)+abs(J)+abs(K);  IJK=[I' J' K'];  % juicy bits
% trace segment from each junction
for i=1:n
    if deg(i)==2 % only start from a junction/terminus
        continue;
    end
    fprintf('  processing node %i of %i\n',i,n);
    junc = [x(i) y(i) z(i)];
    % mark the starting point but remove later
    stack(junc(1),junc(2),junc(3))=10*stack(junc(1),junc(2),junc(3))+1;
    for j=1:deg(i)
        % is the second node already marked? then skip
        now = junc + IJK(neigh{i}(j),:);
        if stack(now(1),now(2),now(3))>10
            continue;
        end
        fprintf('  tracing segment %i from node %i (%i/%i)..\n',ns+1,i,j,deg(i));
        ns=ns+1;
        segs{ns}=nan(300,3);
        segs{ns}(1:2,:)=[junc; now];   nsl=2;
        % special check: 2-node segments
        if stack(now(1),now(2),now(3))~=2
            segs{ns}=segs{ns}(1:2,:);
            continue;
        end
        % mark second stack voxel as having been done
        stack(now(1),now(2),now(3)) = 10*stack(now(1),now(2),now(3));
        % walk to the other end
        go=1;
        while go
            % search 6-neighbour first, then 16, then 26
            vol = stack(now(1)-1:now(1)+1,now(2)-1:now(2)+1,now(3)-1:now(3)+1);
            %vol(2,2,2) = 0; % zero the middle (no need, it should be part of >10)
            vol(vol>10) = 0; % zero the already-visited
            ind = find(vol);
            % now, consider the voxels found in the vicinity
            [yy ci] = min(dist(ind)); % 2-node segment now detected above; ci shouldn't be empty
            if isempty(ci)
                cprintf('red','deg-2 anomaly detected');
                break;
            end
            % see if it's another junction point
            cand=now+IJK(ind(ci),:); 
            if stack(cand(1),cand(2),cand(3))~=2  % stop the seg here
                go = 0;
            else
                stack(cand(1),cand(2),cand(3))=10*stack(cand(1),cand(2),cand(3));
            end
            nsl=nsl+1; segs{ns}(nsl,:)=cand;
            now = cand;
        end
        segs{ns}=segs{ns}(1:nsl,:);
    end
    stack(junc(1),junc(2),junc(3))=(stack(junc(1),junc(2),junc(3))-1)/10;
end
end