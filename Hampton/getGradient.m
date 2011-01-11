function [dzdx, dzdy] = getGradient(ConservativeField)
dzdx = ConservativeField(:,2:end,:) - ConservativeField(:,1:end-1,:);
dzdx = dzdx(2:end,:,:) + dzdx(1:end-1,:,:);
dzdx = dzdx*0.5;

dzdy = ConservativeField(2:end,:,:) - ConservativeField(1:end-1,:,:);
dzdy = dzdy(:,2:end,:) + dzdy(:,1:end-1,:);
dzdy = dzdy*0.5;