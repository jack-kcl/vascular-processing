%% label the vessels as being radial, spiral or canal
%  0: spiral
%  1: radial
%  2: canal
%
%  edit: while we're at it, let's calculate number/volume of maternal canals
function [M segtype mcanals] = label_segments(M, folder)
% see if the files are there
folder = deblank(folder);
fprintf('loading %s...\n',folder);
nrrds=ls([folder filesep '*.nrrd']);
if isempty(nrrds)
    error('Radial and maternal canal masks were not found in %s',folder);
end

% assume the files are named in a standard way
fprintf('reading in radial and canal masks...\n');
canal=nhdr_nrrd_read([folder filesep nrrds(1,:)], 1); % Maternal_Canals.nrrd
radial=nhdr_nrrd_read([folder filesep nrrds(2,:)], 1);  % Radal_ artery.nrrd
canal=single(canal.data);
radial=single(radial.data);
% interpolate the whole lot
rad=interp3(radial,M.nodes(:,2),M.nodes(:,1),M.nodes(:,3));
can=interp3(canal,M.nodes(:,2),M.nodes(:,1),M.nodes(:,3));

segtype = zeros(M.ns,1);
for i=1:M.ns
    % get some nodes from the middle of the seg? or all nodes?
    nds = M.segs{i}.nodes;
    if sum(rad(nds))/length(nds) > 0.5 % why wouldn't it be 0.5, if any of them is true?
        segtype(i)=1; continue;
    end
    if sum(can(nds))/length(nds) > 0.5 % why wouldn't it be 0.5, if any of them is true?
        segtype(i)=2;
    end
end
M=M.setSegmentVariable('segtype',segtype);

%% maternal canals
CC = bwconncomp(canal);
mcanals.num = CC.NumObjects;
mcanals.voxcounts = zeros(mcanals.num,1);
for i=1:mcanals.num
    mcanals.voxcounts(i) = length(CC.PixelIdxList{i});
end
