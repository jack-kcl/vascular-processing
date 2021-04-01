%% given a network with radial/spiral/canal vessels labelled, calculate the
%  global conductance (and local?) from proximal to distal
function [M2 inflow outflow] = calculate_conductance(M2, segtype)

addpath('C:\Users\jl10\Dropbox\CODE\poiseuille');

% separate out the spiral network and find terminal vessels
spsegs = find(segtype==0);
spiral = M2.extractNetwork(segtype==0);
M2 = M2.findJunctions;

% check whether each terminal connects to radial or canal (or none)
spiral = spiral.findBoundary;
bnodes = spiral.getBoundaryNodes('merge');
bn = length(bnodes);

% estimate the proximal/distal status for unknown terminals: 
%     all radial-connected are inlets, all others are outlets  !!!
%inlets = zeros(bn,1); % for each spiral bnode, mark if connected to radial
radials=[];
for i=1:bn
    bnode = spiral.V{1}.data(bnodes(i)); % original node index for M2
    jn = M2.n2j(bnode);
    jdeg = M2.jdeg(jn);
    if jdeg==1; continue; end
    % check if any of the connected segs is a radial artery
    rad = 0;
    for j=1:jdeg
        nd = M2.j2n{jn}(j);
        seg = M2.n2s{nd};
        if segtype(seg)==1
            %rad=1; break;
            radials=[radials; spiral.n2j(bnodes(i))];  % store spiral junc number
            break;
        end
    end
    %inlets(i) = rad;  % will be 1 if connected to radial
end
%radials=find(inlets);

%% build a superstructure (for spiral, not M2) and calculate the conductance
ELS=zeros(spiral.ns,2);
for i=1:spiral.ns
    j1=spiral.segs{i}.nodes(1);
    j2=spiral.segs{i}.nodes(end);
    ELS(i,:) = [spiral.n2j(j1)  spiral.n2j(j2)]; % need to be consecutive unique
    %ELS(i,:) = [j1 j2]; % can't use, nodes are discont at junctions
end
% let's make the ELS consecutive and mark the radial for BC
UELS = unique(ELS);
inds(UELS) = [1:length(UELS)];
CELS = inds(ELS);
BCs = inds(radials);  % global node numbers in poi of inlets
% make the super structure
[spiral LENG] = spiral.calcSegmentLengths;
[spiral RAD] = spiral.averageSegmentVariable('rad');
poi = poiseuille(CELS, RAD, LENG);
%
poi = poi.setViscosity(3.5e-3); % check units
poi = poi.findBoundaryNodes;
for i=1:length(radials)
    %poi = poi.setBC(poi.bnodes(radials(i)), poi.PRES_TYPE, 1.0);
    poi = poi.setBC(BCs(i), poi.PRES_TYPE, 1.0);
end
poi = poi.solve();
poi = poi.calcFlow();
% global conductance is the net flow through the network
poi=poi.calcBoundaryFlow;
% map global inlet indices to bnodes indices
bBCs=zeros(size(BCs));
for j=1:length(BCs)
    bBCs(j) = find(poi.bnodes==BCs(j));
end
inflow = sum(poi.bq(bBCs));
outlets=[1:poi.bnn]; outlets(bBCs)=[];
outflow = sum(poi.bq(outlets));
