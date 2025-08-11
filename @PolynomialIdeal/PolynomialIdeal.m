classdef PolynomialIdeal < handle
    
    properties
        generators
        grobner
    end
    
    methods
        function obj = PolynomialIdeal(generators)
            obj.generators = generators;
            obj.grobner = obj.generators;
        end
        
        function grobnerBasis(obj, reduced)
            if nargin==1
                reduced = 0;
            end
            while obj.extendGrobnerBasis()
            end
            if (reduced)
                % Firstly, make all polynomials of the ideal monic (modify
                % .grobner)
                for i = 1:length(obj.grobner)
                    ci = obj.grobner{i}.leadCoeff();
                    obj.grobner{i} = (1/ci) * obj.grobner{i};
                end
                
                killList = [];
                for i = 1:length(obj.grobner)
                    for j = i+1:length(obj.grobner)
                        lmi = obj.grobner{i}.leadMonomial();
                        lmj = obj.grobner{j}.leadMonomial();
                        if all(lmj <= lmi)
                            killList = [killList; i];
                        end
                    end
                end
                obj.grobner(killList) = [];
            end
        end
        
        function extendedBasis = extendGrobnerBasis(obj)
            n = numel(obj.grobner);
            extendedBasis = 0;
            for i=1:n
                for j=i+1:n
                    sij = obj.grobner{i}.spoly(obj.grobner{j});
                    [~, r] = sij.euclideanDivision(obj.grobner);
                    if ~r.iszero()
                       obj.grobner{end + 1} = r; 
                       extendedBasis = 1;
                       disp(length(obj.grobner))
                    end
                end
            end
        end % -- END extendGrobnerBasis
        
        function flag = ismember(obj, p)
            obj.grobnerBasis();
            [~, r] = p.euclideanDivision(obj.grobner);
            flag = r.iszero();
        end
        
    end % -- END methods
end

