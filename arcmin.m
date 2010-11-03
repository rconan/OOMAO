function radian = arcmin(val)
if ischar(val)
    val = str2double(val);
end
radian = val*constants.arcmin2radian;
end