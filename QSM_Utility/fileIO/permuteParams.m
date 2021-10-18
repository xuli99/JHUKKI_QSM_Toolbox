function Params = permuteParams(Params)
    Params.fov = Params.fov([2,1,3]);
    Params.sizeVol = Params.sizeVol([2,1,3]);
    Params.voxSize = Params.voxSize([2,1,3]);
    if isfield(Params, 'sizeRecAll')
        Params.sizeRecAll = Params.sizeRecAll([2,1,3:length(Params.sizeRecAll)]);
    end
end

