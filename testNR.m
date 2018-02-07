function testNR
	subsize=41;
	numP=6;
	current_folder=pwd;
	% addpath(strcat(current_folder,'\readimxstuff'));
	% [FileName,PathName] = uigetfile('*.im7','Select the images','MultiSelect','on');
	% for i=1:1
		% image_folder=strcat(PathName,FileName(i))
		% image_folder = fullfile( PathName , FileName );
		% I{1}=readimx(image_folder);
	% end
	% F=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
	% [r_g,c_g]=size(F_in);
	% [Xmesh,Ymesh]=meshgrid(1:1:c_g,1:1:r_g);
	% Int = griddedInterpolant(Xmesh',Ymesh',F_in','cubic')
	% % save('Interpolant_func.mat','Int')
	% F_in(1:5,1:5)


	syms dx dy x0 y0 P P1 P2 P3 P4 P5 P6 P7 P8 P9 P10
	A=[(1+P2), P3, P1;
		P5, (P6+1), P4;
		0 0 1]
	X=[dx;dy;1];

	NR3(subsize,numP,A,X)
end