function meshcompare(varargin)
% MESHCOMPARE version 1.3 compares two 2D or 3D arrays, allowing easy visual inspection of
% regions of data (using 'crop'), and other options.
% The figure title provides the line of code to reproduce the first figure's
% surface plot.
% 
% Sometimes it is really hard to compare two arrays due to outliers or large
% offsets. MESHCOMPARE is a straight forward way to compare arrays using
% the matlab MESH/SURF functions (for 3D arrays, a centre slice is taken, see:
% pcslice=0.5; % **
% 
% You can add the string arguments 'crop', 'rotate', and 'subtract' to
% crop, rotate or subtract-the-mean from each field pair before viewing.
% 
% For 3D arrays, you can select 'row','col', or 'z' to specify how the 2D
% mesh is generated. Default is 'z', see direction variable dir=3; % ***
% for example: MESHCOMPARE(aa,aa1,bb,bb1,...,'col') compares aa(:,i,:).
%
% MESHCOMPARE(aa,aa1,bb,bb1,...) compares variables given in input pairs: fields a&a1, b&b1,
% and so on... plotting the first fields (aa,bb) of each pair as a mesh, and the
% second fields (aa1,bb1) as black markers, in separate figures. The titles
% of the figures will be 'aa', and 'bb'.
%
% MESHCOMPARE(aa,aa1,bb,bb1,...,'crop') allows you to select the vertices of
%  a rectangular region to compare in all of the fields pairs (using the
% inbuilt MATLAB function: 'GINPUT').
% MESHCOMPARE(aa,aa1,bb,bb1,...,'rotate') allows you to select a prefered
% viewing orientation for all the field pairs. Very useful to see minute
% differences that are only visible at a certain angle and zoom.
% MESHCOMPARE(aa,aa1,bb,bb1,...,'subtract') subtracts the average values
% of both fields, which is useful if your fields have large scalar offsets
% which make them hard to compare.
%
%
% Example: 2D
%    aa=peaks(50);
%    aa1=aa+randn(size(aa))*0.3; % add some noise
%    bb=real(sqrt(aa));
%    bb1=bb+randn(size(bb))*0.3; % add some noise
%    cc=aa;
%    cc1=real(sin(aa/3)); % compare with similar function
%    meshcomparecrop(aa,aa1, bb,bb1, cc,cc1)
%
% Example: 3D, try the code above replacing line 27 with:
%    aa=flow(50);
%
% Example: crop: try the code above replacing line 33 with:
%    meshcomparecrop(aa,aa1, bb,bb1, cc,cc1,'crop')
%
% Example: crop, subtract, rotate - any combination of these will work:
%    meshcomparecrop(aa,aa1, bb,bb1, cc,cc1,'crop','rotate')
%
% if you have improvements, please email mattmolteno@gmail.com

%% set defaults
% determines percentage through the image that slice is taken (0.5 is 50 %)
pcslice=0.5; % **
% set true to crop the figure
cropfigure=false;
% set true to rotate them to a nice viewpoint
rotatefigure=false;
% subtract mean
subtractmean=false;
% subtract mean
subtractmean=false;
% default: slice rows
dir=3; % ***

mynargin=nargin;

% overwrite defaults if speficied in varargin, then remove from varargin
deletevec=false(1,length(varargin));
for ii=1:length(varargin)
    if ischar(varargin{ii});
        optioni=varargin{ii};
        if isstringthere(optioni,'crop')
            cropfigure=true;
        elseif isstringthere(optioni,'rotate')
            rotatefigure=true;
        elseif isstringthere(optioni,'subtract')
            subtractmean=true;
        elseif isstringthere(optioni,'row')
            dir=1;
        elseif isstringthere(optioni,'col')
            dir=2;
        elseif isstringthere(optioni,'z')
            dir=3;
        end
        deletevec(ii)=true;
    end
end
varargin(deletevec)=[]; % remove from varargin
mynargin=mynargin-sum(deletevec);

% If the remaining inputs is only one, compare with itself:
if mynargin==1
    varargin=[varargin,varargin];
    mynargin=mynargin+1;
end

if rem(mynargin,2)
    disp('Note that MESHCOMPARE has skipped the last variable because it requires an even number of inputs,')
    disp('i.e. input pairs ...aa,aa1,bb,bb1... and so on,');
    disp('see details in description.')
end

% choose example input name:
inputname1=inputname(1);
if isempty(inputname1);inputname1='var1';end

% Set crop area
userinputfig=figure;
if cropfigure
    xfig=varargin{1};
    xfig=midslice_fif(xfig,pcslice,subtractmean,dir);
    imagesc(xfig)
    tempaxis=gca;
    title(['Displaying ''',inputname1,''', select area to look at (vertices of rectangle):'])
    disp(' select area to look at')
    xlabel('your array columns');
    ylabel('rows of your array');
    [cols,rows]=ginput(2);
    cols=sort(round(cols));
    rows=sort(round(rows));
else
    cols=[1;1e6];rows=[1;1e6]; % otherwise select whole images
end
close(userinputfig);

viewfig=figure;
if rotatefigure
    % Similarly: can set the view angle for all fields
    x1=varargin{1};
    x1=midslice_fif(x1,pcslice,subtractmean,dir);
    [colsr,rowsr]=imagebounds_fif(cols,rows,x1);
    x1=x1(rowsr(1):rowsr(2),colsr(1):colsr(2));
    xcomp=varargin{2};
    xcomp=midslice_fif(xcomp,pcslice,subtractmean,dir);
    xcomp=xcomp(rowsr(1):rowsr(2),colsr(1):colsr(2));
    hold on;
    mesh(x1);
    view(3);
    [a,b]=size(xcomp');
    [X,Y]=meshgrid(1:a,1:b);
    plot3(X(:),Y(:),xcomp(:),'k.')
    title([inputname1,' choose view orientation and press enter'])
    disp('choose view orientation and press enter')
    rotate3d on
    pause
    [v1,v2]=view;
    hold off
else
    [v1,v2]=view(3);
end
close(viewfig);

for ii=mynargin-1:-2:1
    
    inputnamei=inputname(ii);
    if isempty(inputnamei);inputnamei=['var',num2str(ii)];end
    inputnameii=inputname(ii+1);
    if isempty(inputnameii);inputnameii=['var',num2str(ii+1)];end
    
    xref=midslice_fif(varargin{ii},pcslice,subtractmean,dir);
    [colsf,rowsf]=imagebounds_fif(cols,rows,xref);
    xref=xref(rowsf(1):rowsf(2),colsf(1):colsf(2));
    xcompare=midslice_fif(varargin{ii+1},pcslice,subtractmean,dir);
    xcompare=xcompare(rowsf(1):rowsf(2),colsf(1):colsf(2));
    if min(size(xref))==1
        error('Field is too small: it might be because the two points you selected in the crop window were too close together.')
    end
    handle_out=figure;
    surf(xref,'edgecolor','interp');
    hold on;
    
    [a,b]=size(xcompare');
    [X,Y]=meshgrid(1:a,1:b);
    plot3(X(:),Y(:),xcompare(:),'k.')
    view(v1,v2)   %uncomment if you want to set view each time (see above)
    alpha(0.8)
%     % % can fiddle with these options if you want things to look fancy
    %     shading flat
    %     shading interp
    %     shading faceted
    %     alpha scaled
    %     view(3)
    axis vis3d square
    
    if ndims(varargin{1})==3
        titleprefix='surf(squeeze(';
        titlesuffix='));';
        switch dir
            case 1
                snum=round(size(varargin{1},1)*pcslice);
                titlestring=['(',num2str(snum),',:,:)'];
                xaxis_text='y (cols)';
                yaxis_text='z (z)';
                zaxis_text='f(y,z)';
            case 2
                snum=round(size(varargin{1},2)*pcslice);
                titlestring=['(:,',num2str(snum),',:)'];
                xaxis_text='x (rows)';
                yaxis_text='z (z)';
                zaxis_text='f(x,z)';
            case 3
                snum=round(size(varargin{1},3)*pcslice);
                titlestring=['(:,:,',num2str(snum),')'];
                xaxis_text='x (cols)';
                yaxis_text='y (rows)';
                zaxis_text='f(x,y)';
        end
    else
        titlestring='';
        titleprefix='surf(';
        titlesuffix=');';
        xaxis_text='x (cols)';
        yaxis_text='y (rows)';
        zaxis_text='f(x,y)';
    end
    title(['>> ',titleprefix,inputnamei,titlestring,titlesuffix])
    xlabel(xaxis_text)
    ylabel(yaxis_text)
    zlabel(zaxis_text)
    try
        legend(inputnamei,inputnameii,'Location','NorthEast')
    catch
        legend(inputnamei,inputnamei,'Location','NorthEast')
        disp('Meshcompare needs even input pairs to compare...')
    end
    hold off

    
end

drawnow


end

function [y]=midslice_fif(x,pcslice,subtractmean,dir)

if ndims(x)==3
    switch dir
        case 1
            y=squeeze(x(round(size(x,1)*pcslice),:,:));
        case 2
            y=squeeze(x(:,round(size(x,2)*pcslice),:));
        case 3
            y=squeeze(x(:,:,round(size(x,3)*pcslice)));
    end
else
    y=x;
end

if subtractmean
    y=y-nanmean(y(:));
end

end

function [cols,rows]=imagebounds_fif(cols,rows,X)
[a,b]=size(X);
cols(cols<1)=1;
rows(rows<1)=1;
cols(cols>b)=b;
rows(rows>a)=a;
end


function [true_false] = isstringthere(string,pattern)
%ISSTRINGTHERE Returns true if pattern is in the ignoring case.

aa=regexp(lower(string),lower(pattern));
if isempty(aa);
    true_false=false; %nothing found
else
    true_false=logical(aa(1)); %found it!
end


end
