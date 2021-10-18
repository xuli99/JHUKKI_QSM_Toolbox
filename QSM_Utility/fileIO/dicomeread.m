function [imag, dicom] = dicomeread(dfile, varargin)
% dicomeread - Read MR enhanced DICOM files with correct order & dimensions
%
%  Uses the Frame Content Sequence, Dimension Index Values to correctly
%  read in MR Enhanced DICOM images with correct n-dimensions.
%  - Adds field MatlabNames to the returned DICOM header containing the
%    labels of the image dimensions. First 2 dimensions are labelled
%    'Column' and 'Row'; other labels come from the DICOM Dimension 
%    Description Label field.
%  - Adds field MatlabIndexes to the returned DICOM header containing
%    the labels of the indexes of each dimension. First 2 dimensions are
%    excluded. Other indexes come from the DICOM field pointed to by the 
%    FunctionalGroupPointer and DimensionIndexPointer. For diffusion scans
%    the indexes are the b-value and directions, for ASL the
%    labeled/unlabelled tag. 
%  - Adds field MatlabItems to the returned DICOM header containing item
%    index for each image. This has the 3rd to end image dimensions. To
%    access frame header info use 
% dicom.PerFrameFunctionalGroupsSequence.
%  - Handles diffusion and other scans with incomplete matricies by
%    combining dimensions E.g. b-value and gradient orientation in case of
%    diffusion, scanning sequence image type in case of B1 maps.
%  - Order property: reorder image dimensions. Arg can be cell string of
%    DICOM dimension names, a row of dimension numbers or 'Prompt' to set
%    interactively.
%  - Combine property: combine dimensions to reduce dimensionality. Arg is
%    a row with 0's and 1's, a dimension with 1 is combined with dimension
%    to the left. Combine is done after re-order.
%  - Noscale Property: don't apply scale slope and intercept to date read.
%
% Example: diffusion scan with multiple dynamics - the order of slices
% in the DICOM file is slices x dynamics x diffusion. The Order arg changes
% this to slices x diffusion x dynamics, then combines slices and diffusion
% to make a 4D file.
% [imag, info] = dicomeread ('diff.dcm', 'Order', [1 2 3 5 4], 'Combine', [0 0 0 1 0]);
%
version = 1.1; %#ok<NASGU>

dorder = 0;
dcombine = 0;
scale = 1;
i = 1;
n = length(varargin);
while i <= n
    arg = varargin{i};
    i = i + 1;
    switch arg
        case 'Combine'
            if i > n
                error ('Dims must be followed by row specifying dimensions to combine');
            end
            dcombine = 1;
            combine = varargin{i};
            i = i + 1;
        case 'Order'
            if i > n
                error ('Order must be followed by row specifying index order');
            end
            dorder = 1;
            order = varargin{i};
            i = i + 1;
        case 'Noscale'
            scale = 0;
        otherwise
            error ('unrecognized property');
    end
end

dicomdict ('set', 'philips-dict.txt');
dicom = dicominfo(dfile);
if ~isfield (dicom, 'PerFrameFunctionalGroupsSequence')
    error 'DICOM file is not MR-Enhanced';
end
imag = dicomread(dicom);
imag = squeeze(imag);
sz = size(imag);
if length(sz) ~= 3
    error ('Image dimensions not correct %d ~= 3', length(sz));
end

% get dimension index seq and fields
ds = dicom.DimensionIndexSequence;
items = fieldnames(ds);
ndims = length(items);
% get per frame functional group and fields
pfgs = dicom.PerFrameFunctionalGroupsSequence;
items = fieldnames(pfgs);
nimag = length(items);
if nimag ~= sz(3)
    error ('Image count not correct %d ~= %d', nimag, sz(3));
end
itemarray = 1:nimag;

% get dimension indexes for all images
dims = zeros(ndims,nimag);
for i = 1:nimag
    dims(:,i) = pfgs.(items{i}).FrameContentSequence.Item_1.DimensionIndexValues;
    ss(:,i) = pfgs.(items{i}).PrivatePerFrameSq.Item_1.MRScaleSlope;
    si(:,i) = pfgs.(items{i}).PrivatePerFrameSq.Item_1.MRScaleIntercept;
end

j = 0;
jj = 0;
comb = [];
for i = 1:ndims
    u = unique(dims(i,:));
    n = length(u);
    % ignore dimensions that have only one value - usually StackID
    if n == 1
        continue
    end
    if u(1) == 1 && u(end) == n
        % dimension goes from 1 to n - no change
        tdims = dims(i,:);
    else
        % else change dim values to be 1 to n
        for k = 1:n
            tdims(dims(i,:) == u(k)) = k;
        end
    end
    
    % if same number of values for each index - no change, store dims,
    % save label, save dim size, save index labels
    t = tabulate(tdims);
    if all(t(:,2) == t(1,2))
        j = j + 1;
        odims(j,:) = tdims;  %#ok<*AGROW>
        nm{j} = replace(ds.(items{i}).DimensionDescriptionLabel, ...
            'Private ', '');
        rdim(j) = n;
        % get labels for each unique index in this dimension
        lbl{j,:} = getLabel(i, ds, pfgs, items, dims, u);
    
    else
        % if not same number of values - this will be combined with one or
        % more other dimensions. save dims, save dim size, save label
        jj = jj + 1;
        comb(jj,:) = tdims;
        nmc{jj} = replace(ds.(items{i}).DimensionDescriptionLabel, ...
            'Private ', '');
        nc(jj) = n;
        % get labels for each unique index in this dimension
        lblc{jj} = getLabel(i, ds, pfgs, items, dims, u);
    
    end
end

% combine dimensions which are not complete into one dimension
if ~isempty(comb)
    % idx is unique index for combined dimensions - may be sparse
    for i = 1:size(comb,2)
        a = num2cell(comb(:,i));
        idx(i) = sub2ind(nc, a{:});
    end
    % get unique values and renumber indexes as 1 to n in odims
    % create labels for each combined index
    u = unique(idx);
    n = length(u);
    j = j + 1;
    for k = 1:n
        odims(j, idx == u(k)) = k;
        lb = '';
        for m = 1:size(comb,1)
            lb = [lb ' ' char(lblc{m}{comb(m,find(idx==u(k),1))})];
        end
        lbs{k} = lb;
    end
    % add to the dim size, save combined label, save index labels in header
    lbl{j,:} = lbs;
    rdim = [rdim n];
    nm{j} = char(join (nmc, '/'));
end

if prod(rdim) ~= nimag
    error('Dimension mismatch %d ~= %d', n, nimag);
end
nm = {'Column', 'Row', nm{:}}; %#ok<CCAT>
n = length(rdim);
oimag = zeros([sz(1:2) rdim]);
oarray = zeros([rdim 1], 'single');
for i = 1:nimag
    if n == 1
        idx = odims(i);
    else
        a = num2cell(odims(:,i));
        idx = sub2ind(rdim, a{:});
    end
    if scale
        oimag(:,:,idx) = (single(imag(:,:,i)) - si(i)) / ss(i);  
    else
        oimag(:,:,idx) = imag(:,:,i);
    end
    oarray(idx) = itemarray(i);
end
imag = oimag;
rdim = size(imag);
itemarray = oarray;

% reorder dimensions if selected
if dorder
    if strcmpi(order, 'prompt')
        if strcmpi(combine, 'prompt')
            [order, combine] = chooseDlg(nm, 1, 1);
        else
            [order, ~] = chooseDlg(nm, 1, 0);
        end
    end
    n = length(order);
    if n == length(rdim)
        err = 0;
        % order is list of names - convert to numbers using dimension
        % labels
        if iscellstr(order) % #ok<ISCLSTR>
            for i = 1:n
                t = find(strcmp(nm,order(i)));
                if isempty(t)
                    fprintf ('Name specified in Order list "%s" not found\nin DICOM Dimension Description Labels\n', order{i});
                    err = 1;
                    break
                end
                ord(i) = t;
            end
        else
            % order is list of dimension numbers
            if ~ (all(order > 0) && all(order <= n) && ...
                    length(unique(order)) == n)
                disp 'Invalid numeric order specified - must be list of numbers representing dimensions for use in permute';
                err = 1;
            else
                ord = order;
                order = nm(ord); % swap labels here - already done for above
            end
        end
        if err
            disp 'Order of dimensions not changed.';
        else
            imag = permute(imag, ord);
            ord = ord(3:end)-2;
            itemarray = permute(itemarray, ord);
            nm = order;
            lbl = lbl(ord);
            rdim = size(imag);
        end
    else
        disp 'Size of order is not same as image dimensions - no permute done'
    end
end

% change number of dimensions if selected
if dcombine
    if strcmpi(combine, 'prompt')
        [~, combine] = chooseDlg(nm, 0, dcombine);
    end
    if length(combine) == length(rdim)
        j = 1;
        for i = 2:length(rdim)
            if combine(i)
                % combine index labels
                for m = 1:rdim(j)
                    for k = 1:rdim(i)
                        lb{m+(k-1)*rdim(j)} = [lbl{j-2}{m} '/' lbl{i-2}{k}];
                    end
                end
                rdim(j) = rdim(j) * rdim(i);
                lbl{j-2} = lb;
                nm{j} = [nm{j} 'x' nm{i}];
            else
                j = j + 1;
                rdim(j) = rdim(i);
                nm(j) = nm(i);
           
            end
        end
        rdim = rdim(1:j);
        nm = nm(1:j);
        imag = reshape (imag, rdim);
        itemarray = reshape (itemarray, rdim(3:end));
    else
        disp 'Size of combine is not same as image dimensions - no combine done'
    end
end
dicom.MatlabNames = nm;
dicom.MatlabIndexes = lbl(1:length(rdim)-2,:);
dicom.MatlabItems = itemarray;
fprintf('Dimension Information: \n');
for dj = 1: length(dicom.MatlabNames)
    fprintf(strcat(nm{dj}, '\n'));
end
end

% get label for all indexes
% i: the dimension to get index labels
% dims: the dimensions index matrix from the DICOM header (ndims x nimages)
% u: the list of unique ordered 1 to n indexes
% ds, pfgs, items: from DICOM header
function lbl = getLabel (i, ds, pfgs, items, dims, u)
% get field names of Dimension Index Pointer and Functional Group
% Pointer for this dimension
d = dicomlookup(dec2hex(ds.(items{i}).DimensionIndexPointer(1)),...
    dec2hex(ds.(items{i}).DimensionIndexPointer(2)));
f = dicomlookup(dec2hex(ds.(items{i}).FunctionalGroupPointer(1)),...
    dec2hex(ds.(items{i}).FunctionalGroupPointer(2)));
for k = 1:length(u)
    % get per-frame functional group item for this image
    pfg = pfgs.(items{find(dims(i,:)==u(k),1)}).(f).Item_1;
    % usual case - return DICOM label field
    if isfield(pfg, d)
        lbl{k} = char(string(pfg.(d)));
    else
        % special handling for diffusion direction field is needed
        if isfield(pfg, 'DiffusionDirectionality')
            % for b=0 or isotropic, there is no direction
            if strcmp(d, 'DiffusionGradientOrientation')
                if strcmp(pfg.DiffusionDirectionality, 'NONE')
                    lbl{k} = '';
                elseif strcmp(pfg.DiffusionDirectionality, 'ISOTROPIC')
                    lbl{k} = 'iso';
                else
                    % for b>0, it is stored in an additional sequence 
                    % inside the diffusion seq as 3 floats
                    lbl{k} = join(string(pfg.DiffusionGradientDirectionSequence.Item_1.DiffusionGradientOrientation), '/');
                end
            end
        else
            % not diffusion and field not found - should never happen
            % another special case may be needed
            lbl{k} = 'unknown';
        end
    end
end
end

function [order, combine] = chooseDlg(nm, dorder, dcombine)
h = length(nm) * 25 + 75;
d = dialog('Position', [300 300 250 h], 'WindowStyle', 'normal',...
    'Name', 'Set Dimensions');
txt{1} = 'Select:';
if dorder
    txt{end+1} = 'Dimension order';
end
if dcombine
    txt{end+1} = 'Dimension(s) to combine with previous';
end
uicontrol('Parent', d, 'Style', 'text',...
    'Position',[20 h-50 210 50], 'String', txt);
for i = 1:length(nm)
    if dorder
        popup(i) = uicontrol('Parent', d, 'Style', 'popup',...
            'Position',[25 h-50-25*i 200 25],...
            'String', nm, 'Value', i, 'UserData', i, ...
            'Callback',@popup_callback);
        if dcombine
            uicontrol('Parent', d, 'Style', 'checkbox',...
                'Position',[225 h-50-25*i 200 25],...
                'String', '', 'Value', 0, 'UserData', i, ...
                'Callback',@cb_callback);
        end
    elseif dcombine
        uicontrol('Parent', d, 'Style', 'checkbox',...
            'Position',[25 h-50-25*i 200 25],...
            'String', nm{i}, 'Value', 0, 'UserData', i, ...
            'Callback',@cb_callback);
    end
end
uicontrol('Parent', d, 'Position', [89 5 70 25],...
    'String',' Close', 'Callback', 'delete(gcf)');
order = nm;
val = 1:length(nm);
combine = zeros (1, length(nm));
% Wait for d to close before running to completion
uiwait(d);

    function popup_callback(cpopup,~)
        popup(val==cpopup.Value).Value = val(cpopup.UserData);
        for j = 1:length(nm)
            val(j) = popup(j).Value;
            order{j} = nm{popup(j).Value};
        end
    end

    function cb_callback(cpopup,~)
        combine(cpopup.UserData) = cpopup.Value;
    end
end