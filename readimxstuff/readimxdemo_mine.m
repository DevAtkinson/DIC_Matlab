[lv_pathstr]= fileparts(which('readimx'));
lv_dirlist = dir([lv_pathstr '/TestImages']);
lv_j=1;
clear A;

       A = readimx([lv_pathstr  '/TestImages/' '3D-PIV one plane.vc7']);
       figure; showimx(A.Frames{1})
       lv_j = lv_j+1;
clear lv_*;