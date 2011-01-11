function [PoissonDist] = randp(lamda,fast,varargin)
%user can decide if they want fast or accurate poisson distributions
error(nargchk(1, 4, nargin))
if nargin == 1
    fast = false;
end
if nargin < 3 || max(size(lamda)) > 1
    [rows col] = size(lamda);
else
    if nargin > 2
        %single amount of photons for whole array
        rows = varargin{1};
        col =  varargin{end}; % allows randp(lamda,rows) to define a square.
    end
end
lamda = round(lamda);
if min(lamda(:)) < 0
    fprintf('can not be negative\n')
    PoissonDist = -1;
    return
end
if max(size(lamda)) == 1
    if lamda > 15
        % 'zero photons unlikely, distribution is same as Gaussian with
        % standard deviation of sqrt(lamda). 0 is 4 sigmas away from lamda when
        % lamda is 16.
        GaussianDist = randn(rows,col);
        GaussianDist = GaussianDist - mean(GaussianDist(:));
        GaussianDist = sqrt(lamda)*GaussianDist/std(GaussianDist(:));
        GaussianDist = GaussianDist + lamda;
        PoissonDist = abs(round(GaussianDist)); %absolute valuse used 'just in case'
    else
        [prob, n] = Poissonian(lamda);
        prob = prob(:);
        n = n(:);
        minProb = max(prob)/10000;
        n(prob < minProb) = [];
        prob(prob < minProb) = [];
        EvenDist = rand(rows,col);
        PoissonDist = zeros(rows,col);
        rHigh = sum(prob(:));
        prob = prob/rHigh;
        rHigh = 1;
        for k = size(n,1):-1:1
            rHigh = rHigh - prob(k);
            PoissonDist(EvenDist > rHigh) = n(k);
            EvenDist(EvenDist > rHigh) = -1;
        end
    end
else
if fast == true
    nP = 1;
else
    nP = 15;
end
    m = [0 1.75 2.91 3.96 4.98 5.99 7:max(lamda(:))]';  %this approach does not work for 1
    s = [0 1.69 1.85 2.06 2.26 2.46 sqrt(7:max(lamda(:)))]';
    %amount of photons varies for array
    [rows col] = size(lamda);
    index = lamda>nP;
    PoissonDist = randn(rows,col);
    EvenDist = rand(rows,col);
 %   GaussianDist(index) = GaussianDist(index) - mean(GaussianDist(index));
    PoissonDist(index) = s(lamda(index)).*PoissonDist(index);
    PoissonDist(index) = PoissonDist(index) + m(lamda(index));
    PoissonDist = abs(round(PoissonDist)); %absolute valuse used 'just in case'
    
    
    
     if isempty(lamda<=nP)
         return
     end
    for p = 1:nP
        % 'zero photons likely, distribution is Poissonian and can not be
        % modelled as a Gaussian.
        temp = zeros(size(lamda(lamda == p)));
        if ~isempty(temp)
            randMap = EvenDist(lamda == p);
            [prob, n] = Poissonian(p);
            prob = prob(:);
            n = n(:);
            minProb = max(prob)/10000;
            n(prob < minProb) = [];
            prob(prob < minProb) = [];
            rHigh = sum(prob(:));
            prob = prob/rHigh;
            rHigh = 1;
            for k = size(n,1):-1:1
                rHigh = rHigh - prob(k);
                temp(randMap > rHigh) = n(k);
                randMap(randMap > rHigh) = -1;
            end
            PoissonDist(lamda == p) = temp;
        end
    end
end

