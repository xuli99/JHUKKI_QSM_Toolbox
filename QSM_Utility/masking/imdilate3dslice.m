function OutMask = imdilate3dslice(InMask, se)

    OutMask = InMask;

    for sliceInd = 1:size(InMask, 3)
        if (~isempty(find(InMask(:,:,sliceInd)>0, 1)))
            OutMask(:,:,sliceInd) = imdilate(InMask(:,:,sliceInd), se);                    
        end

    end