%% take the processed data (segs and R), and clean it up for morphological analysis
%
% name:   name of the mat file to analyse
% plt:    whether to display a matlab plot of nodechains
% export: whether to export a cmgui skelrads exnode file
function M3 = postprocess(name, plt,export)
%function M2 = postprocess(name, plt,export)

if nargin<2
    plt=0;
end
if nargin<3
    export=0;
end

%% load a dataset
load([name '.mat']); %voxelSize=1;   % put um (account for downsizing)

%% plot a network?
% segs: cell array of nodal coordinates
% pt:   junction node coordinates
if plt
    skelrads=zeros(length(R),3); nsk=0;
    figure; leng=zeros(length(segs),1);
    for i=1:length(segs)
        plot3(segs{i}(:,1),segs{i}(:,2),segs{i}(:,3),'b.','markersize',3);
        hold on;
        leng(i) = size(segs{i},1);
        skelrads(nsk+1:nsk+leng(i),:)=segs{i}; nsk=nsk+leng(i);
    end
%     plot3(pt(:,1),pt(:,2),pt(:,3),'r.','markersize',12); axis equal
%     figure; hist(leng,100);
    temp=skelrads; skelrads(:,2)=temp(:,1); skelrads(:,1)=temp(:,2);
end

% export nodes?
if export
    vars{1}.name='radius';
    vars{1}.data=R;
    vars{1}.ncomp=1;
    write_exnode_file([name '_skelrads'],skelrads,nsk,3,vars,1,[1]);
    %write_vtk_file('skelrads',length(pt),0,'tri',1,pt,[], [1],vars);
end

%% damnit, have to put the R back into cell arrays
ns = length(segs);
n=0;
for i=1:ns
    rads{i} = R(n+1:n+size(segs{i},1));
    n = n + size(segs{i},1);
end

%% use meshclass
% we want .nn to be equal (approx..) for every network. Get the total
% length, scale it to 
m = MeshClass.importNodechains(segs, rads);
  m = m.checkSegmentsXi;
  %%m = m.refineNodechain(100, 0.01, 2, segs, rads, [1:2*m.ns]); (n, minLeng, order, nodechain, radchain, segmask)
  %m = m.refineNodechain(50, 0.0, 2, segs, rads, [1:2*m.ns]);  % for segTor2 (scale=leng/50)
  % ---- why is it 2*m.ns?????
  m = m.refineNodeChain2(10000, 2, segs, rads, [1:2*m.ns]); % (N, order, nodechain, radchain, segmask)
nodes=m.nodes;
nodes(:,1)=m.nodes(:,2);
nodes(:,2)=m.nodes(:,1);
m.nodes=nodes; %*voxelSize;

%% clean u by MST
M=m.eliminateRepeatedNodes;
M=M.prepareNetwork;
M=M.findSegments;
M=M.findJunctions;
M=M.calcSegmentLengths;
%[mask M2]=M.MST;  % I've not used the output of MST; it should be run for
%each subnetwork??
M2=M;

%% also remove small unconnected networks
[nsn segsubnet] = find_subnets(M2);
% now consider each subnet and remove if appropriate
keep = ones(M2.ns,1);
leng = M2.getFromSegs('leng');
for i=1:nsn
    ind = find(segsubnet==i);
    % based on length...?
    if sum(leng(ind))<30
        keep(ind)=0;
    end
end
M3 = M2.extractNetwork(keep);
end

%% ---
function [nsn segsubnet] = find_subnets(m)
% just loop over
nsn=0; segsubnet = zeros(m.ns,1);
for i=1:m.ns
    % skip if done
    if segsubnet(i)
        continue; % segment already belongs to a subnet
    end
    % get the second node
    seed = m.segs{i}.nodes(1);
    % get the mask
    mask = m.connectedSubnetwork(seed);
    % update the master mask
    nsn=nsn+1;
    mask=mask*nsn;
    if any(mask & segsubnet)
        error('a segment was found to belong to more than one subnetwork?!');
    end
    segsubnet = segsubnet + mask;  % output contains, per segment, the subnet index
end
end

