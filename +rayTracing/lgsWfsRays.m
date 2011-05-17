src = source;
src.offsetAngle = [zeros(1,8) ; atan2(linspace(-1,1,8)*0.25+12.25,90e3)];
src = src.*abcd('freeSpace',90e3)*abcd('thinLens',90e3)*...
    abcd('freeSpace',100)*abcd('thinLens',4,'offset',12.25)*abcd('freeSpace',5);
figure(101)
clf
rayTrace(src,0,'r')

%%
src1 = source;
src1.offsetAngle = [zeros(1,8) ; atan2(linspace(-1,1,8)*0.25+12.25,85e3)];
src1 = src1.*abcd('freeSpace',85e3)*abcd('thinLens',90e3)*...
    abcd('freeSpace',100)*abcd('thinLens',4,'offset',12.25)*abcd('freeSpace',5);
figure(101)
rayTrace(src1,5e3,'b')

%%
src2 = source;
src2.offsetAngle = [zeros(1,8) ; atan2(linspace(-1,1,8)*0.25+12.25,95e3)];
src2 = src2.*abcd('freeSpace',95e3)*abcd('thinLens',90e3)*...
    abcd('freeSpace',100)*abcd('thinLens',4,'offset',12.25)*abcd('freeSpace',5);
figure(101)
rayTrace(src2,-5e3,'g')