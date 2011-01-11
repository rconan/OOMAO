function phase = PoissonSolve(phase,dPdx,dPdy)
div = getDiv(dPdx,dPdy);
g = .25;
lap = getLaplacian(phase);
phase(2:end-1,2:end-1) = phase(2:end-1,2:end-1) + g*(lap - div);