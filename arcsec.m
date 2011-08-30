function radian = arcsec(val)
if ischar(val)
    val = str2double(val);
end
radian = val*constants.arcsec2radian;
end