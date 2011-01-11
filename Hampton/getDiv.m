function [div curl] = getDiv(dPdx,dPdy)
if max(size(dPdx)) > 1
    [dPdxdx dPdxdy] = getGradient(2*dPdx);
    [dPdydx dPdydy] = getGradient(2*dPdy);
    div = dPdxdx + dPdydy;
    curl = -dPdxdy + dPdydx; %zero for no noise
else
    %not enough data
    div = 0;
    curl = 0;
end
