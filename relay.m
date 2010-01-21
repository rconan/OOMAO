function relay(wave,src)
%% RELAY Source object wave setting
% 
% relay(wave,src) sets the amplitude and phase of the source object based
% on the complex wave input

nSrc = numel(src);
for kSrc = 1:nSrc
    src(kSrc).amplitude = abs(wave);
    src(kSrc).phase     = angle(wave);
end
