%% gather variables to report
function [stats M2] = report_stats_nodes(name,vox)

load([name '.mat']);

% start with the name
stats.name = name;

% grab the inflow
if ~isfield(stats,'inflow')
    stats.inflow = inflow;
end

% number of segments
stats.ns = M2.ns;

% number of subnetworks: find all subnets
M2=M2.findSegments;
M2=M2.buildn2s('fast');
M2=M2.findJunctions;
[stats.nsn stats.segsubnet] = find_subnets(M2);

% (by segment): rad (+-std), length, vol, curvature, tortuosity (types)
% rad
[stats.segRads stats.nodRads] = calc_segrads(M2);  % mean+-sd
% leng
M2 = M2.calcSegmentLengths;
stats.leng = M2.getFromSegs('leng');
M2 = M2.setSegmentVariable('length',stats.leng);
stats.nodLeng = M2.V{M2.getVarIndex('length')}.data;
% vol
[stats.segVol stats.nodVol] = calc_segvol(M2);  % needs checking
% curvature (sum+-sd)
[stats.segCur stats.nodCur] = calc_segcur(M2,vox);

% tortuosity (SOAM, Bullit2008) 
stats.segTor = calc_segtor(M2,vox);

% local DM measure
[stats.segTor2 stats.nodTor2] = calc_segtor2(M2,vox);

% add various fields to M2
M2 = M2.setSegmentVariable('subnet',stats.segsubnet);
M2 = M2.findBoundary;
bnodes = M2.getBoundaryNodes('merge');
bn=zeros(M2.nn,1); bn(bnodes)=5;
M2 = M2.setVariable('terminal',bn); M2.boundary=[];

%% adjust by the effective voxel size (should be in um)
stats.segRads = stats.segRads * vox;
stats.nodRads = stats.nodRads * vox;
stats.leng = stats.leng * vox;
stats.nodLeng = stats.nodLeng * vox;
stats.segVol = stats.segVol * vox^3;
stats.nodVol = stats.nodVol * vox^3;


%fprintf('Note: all metrics are in units of whatever voxel was used!\n');
fprintf('Note: all metrics were adjusted using the effective voxel size %g um\n',vox);

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
    segsubnet = segsubnet + mask;
end
end
%% --
function [segRads nodRads] = calc_segrads(m)
segRads=zeros(m.ns,2); % mean +- sd
ind = m.getVarIndex('rad');
for i=1:m.ns
    r = m.V{ind}.data( m.segs{i}.nodes );
    segRads(i,1) = mean(r);
    segRads(i,2) = std(r);
end
nodRads = m.V{ind}.data;
end
%% --
function [segVol nodVol] = calc_segvol(m)
% estimate segment volumes by truncated cone approximation
segVol = zeros(m.ns,1);
rind = m.getVarIndex('rad');
for i=1:m.ns
    ninds = m.segs{i}.nodes; % ignores the middle node of the element!!
    nodes = m.nodes(ninds,:);
    % element "lengths"
    h = nodes(2:end,:)-nodes(1:end-1,:);
    h = sqrt(h(:,1).^2+h(:,2).^2+h(:,3).^2);
    % volumes
    r = m.V{rind}.data(ninds);
    segVol(i) = sum(pi*(r(1:end-1).^2+r(1:end-1).*r(2:end)+r(2:end).^2).*h)/3;
end
m = m.setSegmentVariable('vol',segVol);
nodVol = m.V{m.getVarIndex('vol')}.data;

end
%% --
function [segCur nodCur] = calc_segcur(m,vox)
m3 = m.prefine(3);
%xi = [0.5];
%xi = [1 2 3 4]/5;
xi = [1:9]/10;
cur = calculate_curvature_1d(m3.ne, m3.nn, m3.nodes*vox, m3.els, m3.B, xi);  % for each el
%cur = rms(cur,2);
cur = max(abs(cur),[],2);
% collect it for a whole segment
segCur = zeros(m.ns,2);
nodCur = zeros(m.ne,1);
for i=1:m.ns
    curs = cur(m.segs{i}.els);
    %segCur(i,:) = [sum(curs) std(curs)];
    segCur(i,:) = [mean(curs) std(curs)];
    % loop for the nodal export (comment out if not needed)
    for j=1:length(m.segs{i}.els)
        nodCur(m.els(m.segs{i}.els(j),:)) = cur(m.segs{i}.els(j));
    end
end
end
%% --
function segTor = calc_segtor(m,vox)
segTor = zeros(m.ns,1);
for i=1:m.ns
    N = length(m.segs{i}.els)*2+1;
    nodes = zeros(N,1);
    nodes(1) = m.segs{i}.els(1,1);
    for j=1:length(m.segs{i}.els)
        % unravel whole segments into nodelist, including middle nodes
        nodes((j-1)*2+1:j*2+1) = m.els(m.segs{i}.els(j),:);
    end
    % apply the formula
    T = m.nodes(nodes(2:end),:)-m.nodes(nodes(1:end-1),:);   T=T *vox;
    IP = zeros(N-3,1); TP = IP; 
    for j=1:length(nodes)-3
        IP(j) = acos(dot(T(j,:)/norm(T(j,:)), T(j+1,:)/norm(T(j+1,:))));
        a = cross(T(j,:),T(j+1,:)); na=norm(a);
        b = cross(T(j+1,:),T(j+2,:)); nb=norm(b);
        if abs(na)<1e-10; na=1; end
        if abs(nb)<1e-10; nb=1; end
        TP(j) = acos(dot( a/na, b/nb ));
    end
    CP = sqrt(IP.^2 + real(TP).^2);
    segTor(i) = sum(CP) / sum(sqrt(T(:,1).^2+T(:,2).^2+T(:,3).^2));
end
segTor = real(segTor);  %% why is this required?!
end
%% --
function [segTor2 nodTor] = calc_segtor2(m,vox)
% use a sliding distance metric (DM) -- apply multiscale?
segTor2 = zeros(m.ns,1);
nodTor = zeros(m.nn,1);
for i=1:m.ns
    % get the full list of nodes
    %ninds = [m.segs{i}.nodes(1):m.segs{i}.nodes(end)]; % Nope. Not safe to assume consecutive numbering
    ninds = m.els(m.segs{i}.els,:)'; temp=ninds(end);
    ninds = ninds(1:2,:); ninds=ninds(:); ninds(end+1)=temp;
    % 
    nn = length(ninds);
    nodes = m.nodes(ninds,:) *vox;
    dist = nodes(2:end,:) - nodes(1:end-1,:);
    dist = sqrt(dist(:,1).^2+dist(:,2).^2+dist(:,3).^2);
    for j=2:nn-1
        % just use steps of +/-1 node? (uncontrolled distance)
        straight = norm(nodes(j+1,:)-nodes(j-1,:));
        nodTor(ninds(j)) = (dist(j-1)+dist(j))/straight;
           nodTor(ninds(j)) = abs(nodTor(ninds(j))-1);  %% hmmm
    end
    segTor2(i) = mean(nodTor(ninds(2:end-1)));
end
end
