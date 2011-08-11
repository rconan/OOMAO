function relay(wave,src)
%% RELAY Source object wave setting
% 
% relay(wave,src) sets the amplitude and phase of the source object based
% on wave input. wave can be either a numerical array of complex amplitudes
% or a cell of the form { amplitude , phase }

nSrc = numel(src);
for kSrc = 1:nSrc
    if iscell(wave)
        src(kSrc).mask      = wave{1}>0;
        src(kSrc).amplitude = wave{1};
        src(kSrc).phase     = wave{2};        
    else
        src(kSrc).mask      = abs(wave)>0;
        src(kSrc).amplitude = abs(wave);
        src(kSrc).phase     = angle(wave);
    end
end
