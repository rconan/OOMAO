classdef gpuLensletArray < lensletArray
% GPULENSLETARRAY Create a lenslet array object
%
% obj = lensletArray(nLenslet) creates a lenslet array object from the
% number of lenslet on one side of the array


    methods

        %% Constructor
        function obj = gpuLensletArray(nLenslet)
            error(nargchk(1, 3, nargin))
            obj = obj@lensletArray(nLenslet);
        end

        function propagateThrough(obj,src)
            
            nSrcStack = size(src,3);
            obj.nArray = size(src,2);
            
            nLensletWavePx = obj.nLensletWavePx;
            nOutWavePx     = 2*nLensletWavePx + rem(nLensletWavePx,2);
            n1             = nLensletWavePx*obj.nLenslet;
            n2             = n1*obj.nArray;
            index = tools.rearrange( [n1,n2] , [nLensletWavePx,nLensletWavePx] );
            
            % phasor to align the spot on the center of the image
            u = (0:(nLensletWavePx-1)).*(1-nLensletWavePx)./nOutWavePx;
            phasor = exp(-1i*pi.*u);
            phasor = phasor.'*phasor;
            
            nLensletTotal = obj.nArray*obj.nLenslet^2;
            w = gpuArray(1:nLensletWavePx);
            index     = gpuArray(index);
            phasor    = gpuArray(single(phasor));
            phasor    = repmat( phasor , [1,1,nLensletTotal] );
%             nOutWavePx     = gpuArray(nOutWavePx);
%             nLensletWavePx = gpuArray(nLensletWavePx);
%             nSrcStack      = gpuArray(nSrcStack);
            intensity = parallel.gpu.GPUArray.zeros(nLensletWavePx,nLensletWavePx,nLensletTotal,'single');
            
            for kSrcStack=1:nSrcStack
                
                m_wave = gpuArray(single(src(1,:,kSrcStack).catWave/nOutWavePx));
                % reshape the lenslet in a 3D array, 1 lenslet in each frame
                m_wave = reshape( m_wave(index), nLensletWavePx , nLensletWavePx , [] );
                m_wave = m_wave.*phasor;
                % and perform the 2D fft
                m_wave( nOutWavePx , nOutWavePx , nLensletTotal ) = 0;
                m_wave = fft2( m_wave);
                
                % apply a field stop the size of the lenslet field of view
                m_wave = m_wave(w,w,:);
                
                % compute the intensity
                intensity = intensity + m_wave.*conj(m_wave);
            end
            clear m_wave
            % reshape intensity into a 2D image
            obj.imagelets        = zeros(n1,n2);
            obj.imagelets(index) = gather(intensity(:));
            obj.imagelets        = obj.imagelets/nSrcStack;%*obj.throughput
            
        end

    end
    
end