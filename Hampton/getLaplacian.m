function L_P = getLaplacian(P)
L_P = conv2(P,[-1 0 -1;0 4 0;-1 0 -1]);
L_P = -L_P(3:end-2,3:end-2);