%% importing & process segment data (ripped from VasEx.m could be a lot cleaner..)
function [knots,nodes,segs] = import_skeleton(folder)
% folder should be the skel folder

%folder='C:\Data\placenta uCT\CFLB3.6i\CFLB3.6i_120417_site 1\down2x\segmentation.seg3dproj\mask-dicom\skel';

readin = fscanf(fopen([folder '\Branch Information.xls']),'%c');

if isempty(str2num(readin(1))) % if first character in the readin is a number then start reading from first character, if not then skip first 77 characters and start readin from character 78
    sc = 78;
else
    sc = 1;
end
fwrite(fopen([folder '\BranchInfo.m'],'wt'),['branchTable = [',readin(sc:size(readin,2)-1),'];']);
fclose('all');

% weird bit
%branchTable = readin(sc:size(readin,2)-1);  % no because it's character type
here = pwd; cd(folder);
BranchInfo; cd(here);

if isempty(str2num(readin(1))) % if first character in the readin is a number then start reading from first character, if not then skip first 77 characters and start readin from character 78
    StartEndTable = branchTable(:,4:9); %start/end coordinates of each branch in skeleton
else
    StartEndTable = branchTable(:,3:8);
end

%Removing same endpoint branches and one/zero voxel length branches
seTable = StartEndTable;
LoopSeg = [];
for row = 1:size(StartEndTable,1)
    if ((seTable(row,1)==seTable(row,4))&&(seTable(row,2)==seTable(row,5))&&(seTable(row,3)==seTable(row,6)))||(branchTable(row,3)<=1)
        LoopSeg = [LoopSeg;row];
    end
end
StartEndTable(LoopSeg,:) = [];

nodesIJ = zeros([size(StartEndTable,1)*2,3],'double');
segsIJ = zeros([size(StartEndTable,1),2],'double');
ni = 1;
for s = 1:size(StartEndTable,1) % iterate over each segment, i.e. each row of the StartEndTable
    segsIJ(s,1:2) = [ni,ni+1];
    nodesIJ(ni,1:3) = StartEndTable(s,1:3); ni = ni+1;
    nodesIJ(ni,1:3) = StartEndTable(s,4:6); ni = ni+1;
end

% read in the skeleton stack
% skelvol = import_png_stack(folder, 'mask-skel');
skelvol = import_png_stack(folder, 'skel');

[knots,nodes,segs] = VesselPathFinder(nodesIJ,segsIJ,skelvol);
