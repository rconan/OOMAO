classdef toeplitzBlockToeplitz < handle
   
    properties
       nRow;
       nCol;
       nBlockRow;
       nBlockCol;
       elements;
       bytes;
       elementsFT;
       lambda;
       mu;
       xi;
       na;
       inverse;
	   tag = 'Toeplitz-Block-Toeplitz';
    end
    
    properties (Access=private)
        log
    end

    methods
        
        function obj = toeplitzBlockToeplitz(nmBlock, blockNM,T)
            obj.nBlockRow = nmBlock(1);
            obj.nBlockCol = nmBlock(2);
            obj.nRow      = blockNM(1);
            obj.nCol      = blockNM(2);
            obj.elements  = T(end:-1:1,end:-1:1);
            
            a = obj.elements;
            a = a(:);
            obj.na = length(a);
            %obj.na = pow2(nextpow2(obj.na));
            obj.elementsFT = fft(a,obj.na);
            
            n = obj.nCol;
            m = obj.nRow;
            [i1,i2] = ndgrid ( 1-m:n-1 ,1-m:n-1 );
            obj.lambda = ( (m + i1 - 1)*(m+n-1) + (m + i2 -1) )';
            [j1,j2] = ndgrid ( 1:n );
            mu_ = (m+n)*(n-1)-(j1-1)*(m+n-1)-(j2-1);
            obj.mu = mu_'+1;
            [j1,j2] = ndgrid ( 1:m );
            xi_ = (m+n)*(m+n-1) - j1*(m+n-1) - j2;
            obj.xi = xi_'+1;

            obj.log = logBook.checkIn(obj);
            %display(obj)
        end
                        
        %% Destructor
        function delete(obj)
            if ~isempty(obj.log)
                checkOut(obj.log,obj)
            end
        end

        function display(obj)
            %% DISPLAY Display object information
            %
            % display(obj) prints information about the toeplitzBlockToeplitz object
            
            fprintf('___ %s ___\n',obj.tag)
            fprintf(' . number of blocks: %dX%d\n',obj.nBlockRow,obj.nBlockCol);
            fprintf(' . size of blocks: %dX%d\n',obj.nRow,obj.nCol);
            fprintf(' . compression factor: %4.0f \n',compressionFactor(obj));
            fprintf('----------------------------------------------------\n')
            
        end        
        
        function t = full(obj,mask)
            T = cell(1,obj.nBlockRow+obj.nBlockCol-1);
            for k=1:length(T)
                T{k} = full( toeplitzMat( obj.nRow , obj.nCol , obj.elements(:,k) ) );
            end
            p = obj.nBlockCol;
            m = obj.nBlockRow;
            x = T;                 % build vector of user data
            cidx = (0:m-1)';
            ridx = p:-1:1;
            t = cidx(:,ones(p,1)) + ridx(ones(m,1),:);  % Toeplitz subscripts
            t = x(t);                                   % actual data    
            t = cell2mat(t);
            if nargin>1 && ~isempty(mask)
                t(mask,:)=[];
                t(:,mask)=[];
            end
        end
        
        function out = transpose(obj)
            % b = a.' computes the non-conjugate transpose of matrix a and
            % returns the result in b.
            T = obj.elements(:,end:-1:1);
            T = T(end:-1:1,:);
            out = toeplitzBlockToeplitz(...
                [obj.nBlockRow,obj.nBlockCol],...
                [obj.nRow,obj.nCol],...
                T);
        end
        
        function out = ctranspose(obj)
            % b = a' computes the complex conjugate transpose of matrix a
            % and returns the result in b
            T = obj.elements(:,end:-1:1);
            T = T(end:-1:1,:);
            T = conj(T);
            out = toeplitzBlockToeplitz(...
                [obj.nBlockRow,obj.nBlockCol],...
                [obj.nRow,obj.nCol],...
                T);
        end
        
        function out = mtimes(obj,b)
            % c = T*b multiplies the matrix by the vector b
            nU = (obj.nBlockRow+obj.nBlockCol-1)*(obj.nRow+obj.nCol-1);
            U = zeros(nU,1);            
            U(obj.mu(:)) = b;
            P = ifft(obj.elementsFT.*fft(U,obj.na));
            out = P(obj.xi(:));
        end
        
        function out = mldivide(obj,c)
            % b = T\c solves Tb = c
            if isempty(obj.inverse)
                add(obj.log,obj,'Computing the inverse...')
                obj.inverse = inv( full(obj) );               
            end
            out = obj.inverse*c;
        end
        
        function P = debugGpu(obj,b)
            obj.na = (obj.nBlockRow+obj.nBlockCol-1)*(obj.nRow+obj.nCol-1);
            U = zeros(obj.na,1);
            U(obj.mu(:)) = b;
            obj.na = pow2(nextpow2(obj.na));
            P = obj.elementsFT.*fft(U,obj.na);
%             P = ifft(obj.elementsFT.*fft(U,na));
%             out = P(obj.xi(:));
        end
        
        function out = maskedMtimes(obj,b,mask)
            out = mask.*mtimes(obj,b);
        end
        
%         function out = size(obj,idx)
%             out = [obj.nBlockRow*obj.nRow,obj.nBlockCol*obj.nCol];
%             if nargin>1
%                 out = out(idx);
%             end
%         end
%         
%         function out = length(obj)
%             out = max( size( obj ) );
%         end
%         
        function out = numel(obj)
            out = 1;
        end
        
        function out = nnz(obj)
            out = (obj.nBlockRow+obj.nBlockCol-1).*(obj.nRow+obj.nCol-1);
        end
        
        function out = compressionFactor(obj)
            out = prod( [obj.nBlockRow*obj.nRow,obj.nBlockCol*obj.nCol] )/nnz(obj);
        end
        
    end
    
end
