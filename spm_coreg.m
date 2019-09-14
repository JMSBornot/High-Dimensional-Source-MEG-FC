function D = spm_coreg(D, val, MEGfids, MRIfids)
persistent S cols

cols = [...
            0         0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840];

S.Din = D;
S.val = val;
S.MEGfids = MEGfids;
S.MRIfids = MRIfids;
S.hedit = zeros(1,12);
S.handles_val = repmat({'00.0'}, [1 12]);
S.mrifids = MRIfids;
S.AZ = [];
S.EL = [];
% only use the fiducial points (ignore digitization points)
S.flag = false;

S.hfig = figure('Name', 'MEG data coregistration', 'NumberTitle', 'off', ...
    'WindowStyle', 'normal', 'Resize', 'off', 'Units', 'normalized', 'Visible', 'off');
set(S.hfig, 'OuterPosition', [0 0 1 1]);

% --- Controls
% texts for row header:
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', 'NAS:', 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.02 0.18 0.03 0.03]);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', 'LPA:', 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.02 0.14 0.03 0.03]);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', 'RPA:', 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.02 0.10 0.03 0.03]);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', 'Head Model Mesh:', 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.02 0.04 0.10 0.02]);
% check boxes:
S.hckbox(1) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'Sagittal', 'Value', 1, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.06 0.186 0.05 0.03], 'Callback', @checkbox_button);
S.hckbox(2) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'Coronal', 'Value', 0, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.11 0.186 0.05 0.03], 'Callback', @checkbox_button);
S.hckbox(3) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'Horizontal', 'Value', 0, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.16 0.186 0.05 0.03], 'Callback', @checkbox_button);
S.hckbox(4) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'Sagittal', 'Value', 0, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.06 0.146 0.05 0.03], 'Callback', @checkbox_button);
S.hckbox(5) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'Coronal', 'Value', 1, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.11 0.146 0.05 0.03], 'Callback', @checkbox_button);
S.hckbox(6) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'Horizontal', 'Value', 1, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.16 0.146 0.05 0.03], 'Callback', @checkbox_button);
S.hckbox(7) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'Sagittal', 'Value', 0, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.06 0.106 0.05 0.03], 'Callback', @checkbox_button);
S.hckbox(8) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'Coronal', 'Value', 0, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.11 0.106 0.05 0.03], 'Callback', @checkbox_button);
S.hckbox(9) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'Horizontal', 'Value', 0, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.16 0.106 0.05 0.03], 'Callback', @checkbox_button);
S.hckbox(10) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'Visible', 'Value', 1, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.12 0.036 0.05 0.03], 'Callback', @checkbox_button);
S.hckbox(11) = uicontrol('Parent', S.hfig, 'Style', 'checkbox', 'String', 'ICP', 'Value', 0, 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.61 0.12 0.05 0.03], 'Callback', @checkbox_button);
% texts for original fiducial MRI coordinate values:
xyz = MRIfids.fid.pnt(strcmp(MRIfids.fid.label,'Nasion'),:);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('X: %2.1f mm', xyz(1)), 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.22 0.19 0.05 0.02]);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('Y: %2.1f mm', xyz(2)), 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.31 0.19 0.05 0.02]);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('Z: %2.1f mm', xyz(3)), 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.40 0.19 0.05 0.02]);
xyz = MRIfids.fid.pnt(strcmp(MRIfids.fid.label,'LPA'),:);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('X: %2.1f mm', xyz(1)), 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.22 0.15 0.05 0.02]);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('Y: %2.1f mm', xyz(2)), 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.31 0.15 0.05 0.02]);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('Z: %2.1f mm', xyz(3)), 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.40 0.15 0.05 0.02]);
xyz = MRIfids.fid.pnt(strcmp(MRIfids.fid.label,'RPA'),:);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('X: %2.1f mm', xyz(1)), 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.22 0.11 0.05 0.02]);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('Y: %2.1f mm', xyz(2)), 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.31 0.11 0.05 0.02]);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('Z: %2.1f mm', xyz(3)), 'FontSize', 12, ...
    'Units','normalized', 'Position', [0.40 0.11 0.05 0.02]);
% texts for modified fiducial MRI coordinate values:
sep = 0.02;
xyz = MRIfids.fid.pnt(strcmp(MRIfids.fid.label,'Nasion'),:);
S.htext(1) = uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('( %2.1f )', xyz(1)), 'FontSize', 12, ...
    'ForegroundColor', cols(7,:), 'Units','normalized', 'Position', [0.22 0.19-sep 0.05 0.02]);
S.htext(2) = uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('( %2.1f )', xyz(2)), 'FontSize', 12, ...
    'ForegroundColor', cols(7,:), 'Units','normalized', 'Position', [0.31 0.19-sep 0.05 0.02]);
S.htext(3) = uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('( %2.1f )', xyz(3)), 'FontSize', 12, ...
    'ForegroundColor', cols(7,:), 'Units','normalized', 'Position', [0.40 0.19-sep 0.05 0.02]);
xyz = MRIfids.fid.pnt(strcmp(MRIfids.fid.label,'LPA'),:);
S.htext(4) = uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('( %2.1f )', xyz(1)), 'FontSize', 12, ...
    'ForegroundColor', cols(7,:), 'Units','normalized', 'Position', [0.22 0.15-sep 0.05 0.02]);
S.htext(5) = uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('( %2.1f )', xyz(2)), 'FontSize', 12, ...
    'ForegroundColor', cols(7,:), 'Units','normalized', 'Position', [0.31 0.15-sep 0.05 0.02]);
S.htext(6) = uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('( %2.1f )', xyz(3)), 'FontSize', 12, ...
    'ForegroundColor', cols(7,:), 'Units','normalized', 'Position', [0.40 0.15-sep 0.05 0.02]);
xyz = MRIfids.fid.pnt(strcmp(MRIfids.fid.label,'RPA'),:);
S.htext(7) = uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('( %2.1f )', xyz(1)), 'FontSize', 12, ...
    'ForegroundColor', cols(7,:), 'Units','normalized', 'Position', [0.22 0.11-sep 0.05 0.02]);
S.htext(8) = uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('( %2.1f )', xyz(2)), 'FontSize', 12, ...
    'ForegroundColor', cols(7,:), 'Units','normalized', 'Position', [0.31 0.11-sep 0.05 0.02]);
S.htext(9) = uicontrol('Parent', S.hfig, 'Style', 'text', 'String', sprintf('( %2.1f )', xyz(3)), 'FontSize', 12, ...
    'ForegroundColor', cols(7,:), 'Units','normalized', 'Position', [0.40 0.11-sep 0.05 0.02]);
% edits for displacements in x,y,z coordinates for modified fiducial points
S.hedit(1) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.27 0.19 0.03 0.02], 'UserData', 1, 'Callback', @edit_button, 'KeyPressFcn', @keyboard_callback);
S.hedit(2) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.36 0.19 0.03 0.02], 'UserData', 2, 'Callback', @edit_button, 'KeyPressFcn', @keyboard_callback);
S.hedit(3) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.45 0.19 0.03 0.02], 'UserData', 3, 'Callback', @edit_button, 'KeyPressFcn', @keyboard_callback);
S.hedit(4) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.27 0.15 0.03 0.02], 'UserData', 4, 'Callback', @edit_button, 'KeyPressFcn', @keyboard_callback);
S.hedit(5) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.36 0.15 0.03 0.02], 'UserData', 5, 'Callback', @edit_button, 'KeyPressFcn', @keyboard_callback);
S.hedit(6) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.45 0.15 0.03 0.02], 'UserData', 6, 'Callback', @edit_button, 'KeyPressFcn', @keyboard_callback);
S.hedit(7) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.27 0.11 0.03 0.02], 'UserData', 7, 'Callback', @edit_button, 'KeyPressFcn', @keyboard_callback);
S.hedit(8) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.36 0.11 0.03 0.02], 'UserData', 8, 'Callback', @edit_button, 'KeyPressFcn', @keyboard_callback);
S.hedit(9) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.45 0.11 0.03 0.02], 'UserData', 9, 'Callback', @edit_button, 'KeyPressFcn', @keyboard_callback);
% controls for rotation input
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', 'Pitch (deg):', 'FontSize', 12, ... % intra-ears axis (roughly pitch or X-axis)
    'Units','normalized', 'Position', [0.22 0.04 0.05 0.02]);
S.hedit(10) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.27 0.04 0.03 0.02], 'UserData', 10, 'Callback', @edit_button);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', 'Roll (deg):', 'FontSize', 12, ... % nose-occipital axis (roughly roll or Y-axis)
    'Units','normalized', 'Position', [0.31 0.04 0.05 0.02]);
S.hedit(11) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.36 0.04 0.03 0.02], 'UserData', 11, 'Callback', @edit_button);
uicontrol('Parent', S.hfig, 'Style', 'text', 'String', 'Yaw (deg):', 'FontSize', 12, ... % ventral-distal axis (roughly yaw or Z-axis)
    'Units','normalized', 'Position', [0.40 0.04 0.05 0.02]);
S.hedit(12) = uicontrol('Parent', S.hfig, 'Style', 'edit', 'String', '00.0', 'FontSize', 10, 'Units', 'normalized', ...
    'Position', [0.45 0.04 0.03 0.02], 'UserData', 12, 'Callback', @edit_button);
% refresh button
uicontrol('Parent', S.hfig, 'Style', 'pushbutton', 'String', 'Refresh', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Units','normalized', 'Position', [0.61 0.05 0.1 0.05], 'Callback', @refresh_button);
% refresh button
uicontrol('Parent', S.hfig, 'Style', 'pushbutton', 'String', 'Reset', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Units','normalized', 'Position', [0.73 0.05 0.1 0.05], 'Callback', @reset_button);
% accept button
uicontrol('Parent', S.hfig, 'Style', 'pushbutton', 'String', 'Accept', 'FontSize', 14, 'FontWeight', 'bold', ...
    'Units','normalized', 'Position', [0.85 0.05 0.1 0.05], 'Callback', @accept_button);

S.hax(1) = axes('Position',[0.01 0.23 0.48 0.76]);
S.hax(2) = axes('Position',[0.51 0.23 0.48 0.76]);
S.ID_SAGITTAL = 1;
S.ID_CORONAL = 2;
S.ID_HORIZONTAL = 3;
S.mesh = spm_eeg_inv_transform_mesh(eye(4), D.inv{val}.mesh);
S.N = nifti(S.mesh.sMRI);
M  = S.N.mat;
d = size(S.N.dat);
% tmp = S.N.mat(1:3,1:3);
% if (sum(sum(abs(tmp - diag(diag(tmp))))) == 0)
%     % matrix physical dimensions are aligned with matrix dimensions
%     S.xyzrange{1} = M(1,1)*(1:d(1)) + M(1,4);
%     S.xyzrange{2} = M(2,2)*(1:d(2)) + M(2,4);
%     S.xyzrange{3} = M(3,3)*(1:d(3)) + M(3,4);
% else


% ndgrid creates grid rationally, i.e. as expected; however interp3
% works only with meshgrid grids. The only difference between ndgrid
% and meshgrid is a permute operation between 1st and 2nd dimension.
[x,y,z] = meshgrid(1:d(1),1:d(2),1:d(3));
tmp = permute(double(S.N.dat), [2 1 3]); % permute data 1st and 2nd dimensions to be compatible with meshgrid format
x1 = M(1,1)*x + M(1,2)*y + M(1,3)*z + M(1,4);
y1 = M(2,1)*x + M(2,2)*y + M(2,3)*z + M(2,4);
z1 = M(3,1)*x + M(3,2)*y + M(3,3)*z + M(3,4);
x1 = x1(:); y1 = y1(:); z1 = z1(:);
xyz = S.mesh.tess_scalp.vert;
S.xyzrange{1} = max(min(xyz(:,1))-5,min(x1)):0.5:min(max(xyz(:,1))+5,max(x1));
if (M(2,2) > 0)
    % this is different from the other to include the nose as much as possible
    S.xyzrange{2} = max(min(xyz(:,2))-15,min(y1)):0.5:min(max(xyz(:,2))+5,max(y1)); 
else
    S.xyzrange{2} = max(min(xyz(:,2))-5,min(y1)):0.5:min(max(xyz(:,2))+15,max(y1));
end
S.xyzrange{3} = max(min(xyz(:,3))-5,min(z1)):0.5:min(max(xyz(:,3))+5,max(z1));
[x1,y1,z1] = ndgrid(S.xyzrange{1},S.xyzrange{2},S.xyzrange{3});
invM = inv(M);
xi = invM(1,1)*x1 + invM(1,2)*y1 + invM(1,3)*z1 + invM(1,4);
yi = invM(2,1)*x1 + invM(2,2)*y1 + invM(2,3)*z1 + invM(2,4);
zi = invM(3,1)*x1 + invM(3,2)*y1 + invM(3,3)*z1 + invM(3,4);
tmp = interp3(x, y, z, tmp, xi(:), yi(:), zi(:));
S.dat = reshape(tmp, size(xi)); 
    
%     error('Unexpected');
% end
[x,y,z,f,xyz,xyznas,vecrot,theta] = deal([]);

% --- Plot figure and make visible this figure
plot_all;
set(S.hfig, 'Visible', 'on');

% --- return output when figure is terminated
uiwait(S.hfig);
% update_xyzFiducials;
if (get(S.hckbox(11),'Value') == 1)
    D = S.Dout;
else
    D = spm_eeg_inv_datareg_ui(D, val, MEGfids, S.mrifids, false);
end
delete(S.hfig);
return

% -------------------------------------------------------------------------------------------------------------------

    % update the slices and mesh plots
    function refresh_button(hObject,eventdata)
        axes(S.hax(1)); cla;
        [S.AZ(1), S.EL(1)] = view;
        axes(S.hax(2)); cla;
        [S.AZ(2), S.EL(2)] = view;        
        plot_all;
    end

    % reset the slices and mesh plots to its initial position and
    % orientation
    function reset_button(hObject,eventdata)
        S.mrifids = S.MRIfids;
        S.AZ = [];
        S.EL = [];
        set(S.hedit,'String','00.0');
        tmp = S.mrifids.fid.pnt'; tmp = tmp(:);
        for it = 1:9
            set(S.htext(it), 'String', sprintf('( %2.1f )', tmp(it)));
        end
        S.handles_val = repmat({'00.0'}, [1 12]);
        set(S.hckbox([1 5 6 10]), 'Value', 1);
        set(S.hckbox([2:4 7:9 11]), 'Value', 0);
        plot_all;
    end

    % accept the changes and close the figure
    function accept_button(hObject,eventdata) %#ok<*INUSD>
        if strcmp(S.lastaction,'updated')
            uiresume(S.hfig);
        else
            errordlg('Please click first on Refresh button to check coregistration.', 'Error Dialog', 'modal');
        end
    end

    % callback when subject interact with checkbox
    function checkbox_button(hObject,eventdata)
        axes(S.hax(1)); cla; % colormap gray
        [S.AZ(1), S.EL(1)] = view;
        plot_slices;
        if strcmp(get(hObject,'String'), 'ICP')
            S.lastaction = 'recompute';
        end
    end

    % callback when subject interact with edit buttons
    function edit_button(hObject,eventdata)
        num = str2num(get(hObject, 'String'));
        if isempty(num)
            errordlg('Please enter a real number.', 'Error Dialog', 'modal');
            set(hObject, 'String', '00.0');
        end
        id = get(hObject, 'UserData');
        switch id
            case {1 2 3 4 5 6 7 8 9}
                if (num < -99.9) || (num > 99.9)
                    errordlg('Allowed displacement must be between -99.9 and 99.9 millimiters.', 'Error Dialog', 'modal');
                    set(hObject, 'String', S.handles_val{id});
                else
                    S.handles_val{id} = get(hObject, 'String');
                end
            case {10 11 12}
                if (num <= -180) || (num >= 180)
                    errordlg('Please insert an angle in the main sector: -180 to 180 degrees.', 'Error Dialog', 'modal');
                    set(hObject, 'String', S.handles_val{id});
                else
                    S.handles_val{id} = get(hObject, 'String');
                end
        end
        S.lastaction = 'recompute';
    end

    % how to manage keyboard inputs when control is in edits for
    % translation input?
    function keyboard_callback(hObject, eventdata)
        set(S.hckbox(11),'Value',0);
        id = get(hObject, 'UserData');
        num = str2num(S.handles_val{id});
        switch eventdata.Key
            case 'leftarrow'
                num = num - 1;
            case 'rightarrow'
                num = num + 1;
        end
        if (num < -99.9) || (num > 99.9)
            errordlg('Allowed displacement must be between -99.9 and 99.9 millimiters.', 'Error Dialog', 'modal');
        else
            S.handles_val{id} = sprintf('%2.1f', num);
            set(hObject, 'String', S.handles_val{id});
            refresh_button;
        end
    end

    % update fiducial coordinates
    function update_xyzFiducials
        S.lastaction = 'updated';
        S.mrifids = S.MRIfids;
        xyznas = S.mrifids.fid.pnt(strcmp(S.mrifids.fid.label,'Nasion'),:);
        theta = str2num(S.handles_val{10}); % rotation around the "X" axis
        if (theta ~= 0)
            % first, rotation around (pitch) axis?
            theta = pi*theta/180; % convert from degrees to radians
            vecrot = S.mrifids.fid.pnt(strcmp(S.mrifids.fid.label,'LPA'),:) - S.mrifids.fid.pnt(strcmp(S.mrifids.fid.label,'RPA'),:);
            update_xyzFiducials_rotation('LPA');
            update_xyzFiducials_rotation('RPA');
        end
        theta = str2num(S.handles_val{11}); % rotation around the "Y" axis
        if (theta ~= 0)
            % second, rotation around (roll) axis?
            theta = pi*theta/180; % convert from degrees to radians
            vecrot = xyznas - 0.5*(S.mrifids.fid.pnt(strcmp(S.mrifids.fid.label,'LPA'),:)+S.mrifids.fid.pnt(strcmp(S.mrifids.fid.label,'RPA'),:));
            update_xyzFiducials_rotation('LPA');
            update_xyzFiducials_rotation('RPA');
        end
        theta = str2num(S.handles_val{12}); % rotation around the "Z" axis
        if (theta ~= 0)
            % third, rotation around (yaw) axis?
            theta = pi*theta/180; % convert from degrees to radians
            vecrot = cross(S.mrifids.fid.pnt(strcmp(S.mrifids.fid.label,'LPA'),:)-xyznas, S.mrifids.fid.pnt(strcmp(S.mrifids.fid.label,'RPA'),:)-xyznas);
            update_xyzFiducials_rotation('LPA');
            update_xyzFiducials_rotation('RPA');
        end
        % finally, apply translation
        update_xyzFiducials_translation;
        % fill up the text fields with updated coordinates
        tmp = S.mrifids.fid.pnt'; tmp = tmp(:);
        for it = 1:9
            set(S.htext(it), 'String', sprintf('( %2.1f )', tmp(it)));
        end
    end

    % update translation coordinates accordingly with rotation
    function update_xyzFiducials_rotation(strfid)
        flag = strcmp(S.mrifids.fid.label,strfid);
        vec = S.mrifids.fid.pnt(flag,:) - xyznas;
        vec = axis_rot(vec, vecrot, theta);
        S.mrifids.fid.pnt(flag,:) = xyznas' + vec;
    end

    % apply translation coordinates for fiducial points
    function update_xyzFiducials_translation
        xyz = [str2num(get(S.hedit(1),'String')) str2num(get(S.hedit(2),'String')) str2num(get(S.hedit(3),'String'))]; %#ok<*ST2NM>
        flag = strcmp(S.mrifids.fid.label, 'Nasion');
        S.mrifids.fid.pnt(flag,:) = S.mrifids.fid.pnt(flag,:) + xyz;
        xyz = [str2num(get(S.hedit(4),'String')) str2num(get(S.hedit(5),'String')) str2num(get(S.hedit(6),'String'))];
        flag = strcmp(S.mrifids.fid.label, 'LPA');
        S.mrifids.fid.pnt(flag,:) = S.mrifids.fid.pnt(flag,:) + xyz;
        xyz = [str2num(get(S.hedit(7),'String')) str2num(get(S.hedit(8),'String')) str2num(get(S.hedit(9),'String'))];
        flag = strcmp(S.mrifids.fid.label, 'RPA');
        S.mrifids.fid.pnt(flag,:) = S.mrifids.fid.pnt(flag,:) + xyz;
    end

    % prepare coordinate and function values to plot slices with surf
    function prepare_slices(id)
        % correspondingly either id=1 (sagittal), or id=2 (coronal), or
        % id=3 (horizontal)
        [~,ind] = min(abs(S.xyzrange{id} - xyz(id)));
        switch id
            case 1
                f = squeeze(S.dat(ind,:,:)); % sagittal
                [x,y,z] = ndgrid(S.xyzrange{id}(ind), S.xyzrange{2}, S.xyzrange{3});
            case 2
                f = squeeze(S.dat(:,ind,:)); % coronal
                [x,y,z] = ndgrid(S.xyzrange{1}, S.xyzrange{id}(ind), S.xyzrange{3});
            case 3
                f = S.dat(:,:,ind); % horizontal
                [x,y,z] = ndgrid(S.xyzrange{1}, S.xyzrange{2}, S.xyzrange{id}(ind));
        end
        x = squeeze(x);
        y = squeeze(y);
        z = squeeze(z);
    end

    % plot all the outputs
    function plot_all
        update_xyzFiducials;
        if (get(S.hckbox(11),'Value') == 1)
            Finter = spm('FnUIsetup','MEEG/MRI coregistration', 0);
            Fgraph  = spm_figure('GetWin','Graphics');
            spm_figure('Clear',Fgraph);
            % first ICP registration
            S.Dout = spm_eeg_inv_datareg_ui(S.Din, S.val, MEGfids, S.mrifids, true);
            datareg = S.Dout.inv{1}.datareg;
            flag = strcmp({datareg.modality}, 'MEG');
            datareg = datareg(flag);
            mrifids = datareg.fid_mri;
            % then rigid registration
            cfg = [];
            cfg.sourcefid = mrifids;
            cfg.targetfid = S.mrifids;
            cfg.useheadshape = false;
            S.M1 = spm_eeg_inv_datareg(cfg);
            close(Finter);
            close(Fgraph);
        else
            cfg = [];
            cfg.sourcefid = S.MEGfids;
            cfg.targetfid = S.mrifids;
            cfg.useheadshape = false;
            S.M1 = spm_eeg_inv_datareg(cfg); % rigid registration
        end
        axes(S.hax(1)); cla;
        plot_slices;
        axes(S.hax(2)); cla;
        plot_coreg;
    end

    % plot slices + fiducial points
    function plot_slices
        meegfid = ft_transform_headshape(S.M1, S.MEGfids);
        hold on;
        % localise the slices for: Nasion
        xyz = S.mrifids.fid.pnt(strcmp(S.mrifids.fid.label,'Nasion'),:);
        plot3(xyz(1), xyz(2), xyz(3), 'd', 'MarkerFaceColor', cols(4,:), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
        if (get(S.hckbox(1),'Value') == 1)
            prepare_slices(S.ID_SAGITTAL);
            surf(x, y, z, f, 'EdgeColor', 'none'); shading interp;
        end
        if (get(S.hckbox(2),'Value') == 1)
            prepare_slices(S.ID_CORONAL);
            surf(x, y, z, f, 'EdgeColor', 'none'); shading interp;
        end
        if (get(S.hckbox(3),'Value') == 1)
            prepare_slices(S.ID_HORIZONTAL);
            surf(x, y, z, f, 'EdgeColor', 'none'); shading interp;
        end
        % LPA
        xyz = S.mrifids.fid.pnt(strcmp(S.mrifids.fid.label,'LPA'),:);
        plot3(xyz(1), xyz(2), xyz(3), 'd', 'MarkerFaceColor', cols(4,:), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
        if (get(S.hckbox(4),'Value') == 1)
            prepare_slices(S.ID_SAGITTAL);
            surf(x, y, z, f, 'EdgeColor', 'none'); shading interp;
        end
        if (get(S.hckbox(5),'Value') == 1)
            prepare_slices(S.ID_CORONAL);
            surf(x, y, z, f, 'EdgeColor', 'none'); shading interp;
        end
        if (get(S.hckbox(6),'Value') == 1)
            prepare_slices(S.ID_HORIZONTAL);
            surf(x, y, z, f, 'EdgeColor', 'none'); shading interp;
        end
        % RPA
        xyz = S.mrifids.fid.pnt(strcmp(S.mrifids.fid.label,'RPA'),:);
        plot3(xyz(1), xyz(2), xyz(3), 'd', 'MarkerFaceColor', cols(4,:), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
        if (get(S.hckbox(7),'Value') == 1)
            prepare_slices(S.ID_SAGITTAL);
            surf(x, y, z, f, 'EdgeColor', 'none'); shading interp;
        end
        if (get(S.hckbox(8),'Value') == 1)
            prepare_slices(S.ID_CORONAL);
            surf(x, y, z, f, 'EdgeColor', 'none'); shading interp;
        end
        if (get(S.hckbox(9),'Value') == 1)
            prepare_slices(S.ID_HORIZONTAL);
            surf(x, y, z, f, 'EdgeColor', 'none'); shading interp;
        end
        if (get(S.hckbox(10),'Value') == 1)
            patch('vertices', S.mesh.tess_scalp.vert, 'faces', S.mesh.tess_scalp.face, ...
                'EdgeColor', cols(3,:), 'EdgeAlpha', 0.3, 'FaceColor', 'none');
        end
        plot3(meegfid.pnt(:,1), meegfid.pnt(:,2), meegfid.pnt(:,3), 'o', ...
            'MarkerFaceColor', cols(6,:), 'MarkerSize', 4, 'MarkerEdgeColor', 'k');
        plot3(meegfid.fid.pnt(:,1), meegfid.fid.pnt(:,2), meegfid.fid.pnt(:,3), 'o', ...
            'MarkerFaceColor', cols(6,:), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
        hold off;
        rotate3d on; axis image off;
        if isempty(S.AZ)
            view(-124, 10);
        else
            view(S.AZ(1), S.EL(1));
        end
    end

    % plot co-registration: brain meshes + fiducial points
    function plot_coreg
        mrifid = S.mrifids;
        meegfid = ft_transform_headshape(S.M1, S.MEGfids);
        mesh = S.mesh;
        grad = ft_transform_headshape(S.M1, S.Din.sensors('MEG'));
        hold on
        % plot
        ft_plot_sens(grad);
        patch('vertices', mesh.tess_ctx.vert, 'faces', mesh.tess_ctx.face, ...
            'EdgeColor', cols(1,:), 'FaceColor', cols(1,:));
        patch('vertices', mesh.tess_iskull.vert, 'faces', mesh.tess_iskull.face, ...
            'EdgeColor', cols(2,:), 'FaceColor', 'none');
        patch('vertices', mesh.tess_scalp.vert, 'faces', mesh.tess_scalp.face, ...
            'EdgeColor', cols(3,:), 'FaceColor', 'none');
        plot3(meegfid.pnt(:,1), meegfid.pnt(:,2), meegfid.pnt(:,3), 'o', ...
            'MarkerFaceColor', cols(6,:), 'MarkerSize', 4, 'MarkerEdgeColor', 'k');
        plot3(meegfid.fid.pnt(:,1), meegfid.fid.pnt(:,2), meegfid.fid.pnt(:,3), 'o', ...
            'MarkerFaceColor', cols(6,:), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
        plot3(mrifid.fid.pnt(:,1), mrifid.fid.pnt(:,2), mrifid.fid.pnt(:,3), 'd', ...
            'MarkerFaceColor', cols(4,:), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
        hold off
        rotate3d on;
        axis image;
        if isempty(S.AZ)
            view(-124, 10);
        else
            view(S.AZ(2), S.EL(2));
        end
    end
end