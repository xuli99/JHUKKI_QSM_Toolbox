function tf = isDICOM(filename)

% Open the file -- You can do this in the native endian type.
fid = fopen(filename, 'r');

fseek(fid, 128, 'bof');

if (isequal(fread(fid, 4, 'char=>char')', 'DICM'))
    
    % It has the form of a compliant DICOM file.
    tf = true;
    
else
    
    % It may be a DICOM file without the standard header.
    fseek(fid, 0, 'bof');
    
    tag = fread(fid, 2, 'uint32')';
    
    if ((isequal(tag, [8 0]) || isequal(tag, [134217728 0])) || ...
            (isequal(tag, [8 4]) || isequal(tag, [134217728 67108864])) || isequal(tag, [327688 10]))
        
        % The first eight bytes look like a typical first tag.
        tf = true;
        
    else
        
        % It could be a DICOM file, but it's hard to say.
        tf = false;
        
    end
    
end

fclose(fid);