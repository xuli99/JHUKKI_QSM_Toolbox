function [finalmask] = cjlabel( maskin,groups )

finalmask = zeros([size(maskin) groups]);

for groupnum = 1:groups
    bmask = zeros([size(maskin)]);

    for ii =1:size(maskin,3)
        L = bwlabel( maskin(:,:,ii) );
        if L == 0
            bmask(:,:,ii) = maskin(:,:,ii);
        else
            [counts, bins] = hist( L(:), 0:max(L(:)));
            [junk, inds] = sort( counts );
            if (length(inds) - groupnum > 0)
                bmask(:,:,ii) =  L == bins(inds(end-groupnum));
            end
        end
    end
    finalmask(:,:,:,groupnum) = bmask;

end



