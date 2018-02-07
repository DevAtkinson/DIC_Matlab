function mask=createMask(F_in)
	scrsz = get(0,'ScreenSize');
	width_of_gui=950;
	height_of_gui=600;
	handles.fig = figure('MenuBar','None','Position',[(scrsz(3)-width_of_gui)/2 (scrsz(4)-height_of_gui)/2 width_of_gui height_of_gui]);
	axes1 = axes('Parent',handles.fig ,'Layer' ,'Top','units','pixels','Position',[300 100 400 400]);
	imagesc(F_in)
	button_width=80;
    button_height=25;
    roicount=0;
    roi(1).handle = [];
    roi(1).act = [];
    roi(1).coord = [];
    roi(1).shape = [];
    roi(1).clude = 'include';
    polygon_button=uicontrol('Parent',handles.fig ,'style','pushbutton','string','Polygon','Position',[5 15 button_width button_height],'Callback',@Poly);
    ellipse_button=uicontrol('Parent',handles.fig ,'style','pushbutton','string','Ellipse','Position',[5 45 button_width button_height],'Callback',@Ellip);
    rect_button=uicontrol('Parent',handles.fig ,'style','pushbutton','string','Rectangle','Position',[5 75 button_width button_height],'Callback',@Rect);
    fin=uicontrol('Parent',handles.fig ,'style','pushbutton','string','Done','Position',[5 400 button_width button_height],'Callback',@finished);
	uiwait

	function Poly(hObject, eventdata, handles)
        % hObject    handle to MHData_Select (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        % update handles
        handles = guidata(hObject);
        % set the correct axes to apply the mask to
        axes(axes1);
        fprintf('~ Creating polygon mask\n');
        %create dragable rectangle
        roicount = roicount+1;
        roi(roicount).handle = impoly;
        set( get(roi(roicount).handle,'Children'),....
                'UIContextMenu','' );
        roi(roicount).act='mask';
        roi(roicount).shape = 1;
        roi(roicount).coord = getPosition(roi(roicount).handle);
        roi(roicount).clude = 'include';
        
        setColor(roi(roicount).handle,'g');
        addNewPositionCallback(roi(roicount).handle,@(pos)callbackroi(roicount));
        handles.roi(roicount).coord=roi(roicount).coord;
        handles.roi=roi;

        guidata(hObject, handles);
        uiwait
    end

    function Rect(hObject, eventdata, handles)
        % hObject    handle to MHData_Select (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        % update handles
        handles = guidata(hObject);
        % set the correct axes to apply the mask to
        axes(axes1);
        fprintf('~ Creating rectangular mask\n');
        %create dragable rectangle
        roicount = roicount+1;
        roi(roicount).handle = imrect;
        set( get(roi(roicount).handle,'Children'),....
                'UIContextMenu','' );
        roi(roicount).act='mask';
        roi(roicount).shape = 2;
        roi(roicount).coord = getPosition(roi(roicount).handle);
        roi(roicount).clude = 'include';
        
        setColor(roi(roicount).handle,'g');
        addNewPositionCallback(roi(roicount).handle,@(pos)callbackroi(roicount));
        handles.roi(roicount).coord=roi(roicount).coord;
        handles.roi=roi;

        guidata(hObject, handles);
        uiwait
    end

    function Ellip(hObject, eventdata, handles)
        % hObject    handle to MHData_Select (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        % update handles
        handles = guidata(hObject);
        % set the correct axes to apply the mask to
        axes(axes1);
        fprintf('~ Creating rectangular mask\n');
        %create dragable rectangle
        roicount = roicount+1;
        roi(roicount).handle = imellipse;
        set( get(roi(roicount).handle,'Children'),....
                'UIContextMenu','' );
        roi(roicount).act='mask';
        roi(roicount).shape = 3;
        roi(roicount).coord = getPosition(roi(roicount).handle);
        roi(roicount).clude = 'include';
        
        setColor(roi(roicount).handle,'g');
        addNewPositionCallback(roi(roicount).handle,@(pos)callbackroi(roicount));
        handles.roi(roicount).coord=roi(roicount).coord;
        handles.roi=roi;

        guidata(hObject, handles);
        uiwait
    end

    function finished(hObject,eventdata,handles)
    	handles = guidata(hObject);
    	mask=ones(size(F_in));
    	for i=1:roicount
    		if handles.roi(i).shape==1
    			imagesc(F_in)
    			it=impoly(axes1,handles.roi(i).coord);
    			mask_temp=createMask(it);
    		elseif handles.roi(i).shape==2
    			imagesc(F_in)
    			it=imrect(axes1,handles.roi(i).coord);
    			mask_temp=createMask(it);
    		elseif handles.roi(i).shape==3
    			imagesc(F_in)
    			it=imellipse(axes1,handles.roi(i).coord);
    			mask_temp=createMask(it);
    		end
   %  		if handles.roi(i).clude=='include'
			% 	mask_temp=double(mask_temp);
			% 	mask_temp(mask_temp==0)=NaN;
			% end
			if handles.roi(i).clude=='exclude'
				mask_temp=not(mask_temp);
				% mask_temp=double(mask_temp);
				% mask_temp(mask_temp==0)=NaN;
			end
    		mask=mask.*mask_temp;
    	end
    	guidata(hObject, handles);
    	uiresume
		close
    end

    function callbackroi( roicount )
        % CALLBACKROI resets properties of current selected roi.
        %   CALLBACKROI(roicount) resets position, mask, colour and label for
        %   selected roi.
        
        % Set new positions
        roi(roicount).coord = getPosition(roi(roicount).handle);

        % Set new UI menu (right click).
        rightclickhandle = uicontextmenu('parent',handles.fig);
        uimenu(rightclickhandle,'Label','Include','Callback',@rightclickmenu);
        uimenu(rightclickhandle,'Label','Exclude','Callback',@rightclickmenu);
        uimenu(rightclickhandle,'Label','Delete' ,'Callback',@rightclickmenu);
        uimenu(rightclickhandle,'Label','Mask range' ,'Callback',@rightclickmenu);
        set(get(roi(roicount).handle,'Children'),'UIContextMenu',rightclickhandle);
        
        function rightclickmenu( source,~,~ ) %create a rightclickmenu with include/exclude/mask range options
            % RIGHTCLICKMENU options called in uicontextmenu and uimenu.
            %   RIGHTCLICKMENU(source,...,handle) are Matlab defined varargins.
            
            switch source.Label
                case 'Include'
                    setColor(roi(roicount).handle,'g')
                    roi(roicount).clude = 'include';
                case 'Exclude'
                    setColor(roi(roicount).handle,'r')
                    roi(roicount).clude = 'exclude';
                case 'Delete'
                    delete(roi(roicount).handle);
                    roi(roicount).act ='removed';
                case 'Mask range'
                    if isDVC
                        if (roi(roicount).Xaxis=='PosX')&(roi(roicount).Yaxis=='PosY')
                            [mask_handle,out1,out2]=mask_range_gui(zmax)
                            max_check=zmax;
                        elseif (roi(roicount).Xaxis=='PosX')&(roi(roicount).Yaxis=='PosZ')
                            [mask_handle,out1,out2]=mask_range_gui(ymax)
                            max_check=ymax;
                        elseif (roi(roicount).Xaxis=='PosY')&(roi(roicount).Yaxis=='PosZ')
                            [mask_handle,out1,out2]=mask_range_gui(xmax)
                            max_check=xmax;
                        end
                        val1=out1;
                        val2=out2;
                        if (val1<=max_check)&(val2<=max_check)&(val1>0)&(val2>0)
                            roi(roicount).mask_range=sort([val1 val2]);
                        else
                            fprintf('~ Error: The selected range does not fall within the range of the data\n');
                        end
                    else 
                        fprintf('Applying a mask range is only available for DVC and FEM data\n');
                    end
            end
        end
    end
end