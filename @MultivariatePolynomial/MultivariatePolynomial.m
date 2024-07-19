classdef MultivariatePolynomial < handle
    %MULTIVARIATEPOLYNOMIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        matrixData
        ord
        varNames
    end
    
    methods
        function obj = MultivariatePolynomial(matrixData, ord, varNames)
            %MULTIVARIATEPOLYNOMIAL Construct an instance of this class
            %   Detailed explanation goes here
            if nargin <= 1
                ord = @(s1, s2) lex(s1, s2);
            end
            obj.ord = ord;
            obj.matrixData = monomial_sort(matrixData, ord);
            obj.collectMonomials()
            if nargin <= 2
                obj.varNames = cell(obj.numMonomials, 1);
                for i=1:obj.numMonomials
                    obj.varNames{i} = ['x', num2str(i)];
                end
            else
                obj.varNames = varNames;
            end
        end
        
        function n = numMonomials(obj)
            n = size(obj.matrixData, 1);
        end
        
        function n = numIndeterminates(obj)
            n = size(obj.matrixData, 2) - 1;
        end
        
        function lm = leadMonomial(obj)
            lm = obj.matrixData(1, 1:end-1);
        end
        
        function lc = leadCoeff(obj)
            lc = obj.matrixData(1, end);
        end
        
        function lt = leadTerm(obj)
            lt = obj.matrixData(1, :);
        end
        
        function f = iszzero(obj)
            f = isempty(obj.matrixData) || all(obj.matrixData(:, end)==0);
        end
        
        function f = ismonomial(obj)
            f = (obj.numMonomials == 1);
        end
        
        function r = divide_lead_term_by_lead_term_of(obj, other)
            div = monomial_divide(obj.leadMonomial, other.leadMonomial);
            if isempty(div)
                r = [];
                return;
            end
            r = MultivariatePolynomial([monomial_divide(obj.leadMonomial, other.leadMonomial) obj.leadCoeff/other.leadCoeff], obj.ord, obj.varNames);
        end
        
        function collectMonomials(obj)
            nM = obj.numMonomials();
            iota_cache = zeros(nM, 1);
            iota_cache(1) = 1;
            for i=1:nM-1
                if all(obj.matrixData(i, 1:end-1) == obj.matrixData(i+1, 1:end-1))
                    iota_cache(i+1, 1) = iota_cache(i, 1);
                    obj.matrixData(i+1, end) = ...
                        obj.matrixData(i+1, end) + obj.matrixData(i, end); % accummulate
                else
                    iota_cache(i+1, 1) = i+1;
                end
            end
            taboo = flipud(diff(flipud(iota_cache)));
            obj.matrixData(taboo==0, :) = [];
            obj.matrixData(obj.matrixData(:, end)==0, :) = [];
        end % END collectMonomials
        
        function r = plus(obj1,obj2)
            M1 = obj1.matrixData;
            M2 = obj2.matrixData;
            sumMat = [M1;M2];
            r = MultivariatePolynomial(sumMat, obj1.ord, obj1.varNames);
        end % END plus
        
        function r = minus(obj1, obj2)
            M1 = obj1.matrixData;
            M2 = obj2.matrixData;
            M2(:, end) = -M2(:, end);
            sumMat = [M1;M2];
            r = MultivariatePolynomial(sumMat, obj1.ord, obj1.varNames);
        end % END plus
        
        function r = mtimes(obj1,obj2)
            %Polynomial-polynomial and scalar-polynomial multiplication
            if isa(obj1,'double') && isa(obj2,'MultivariatePolynomial')
                t = zeros(1, obj2.numIndeterminates()+1);
                t(end) = obj1;
                obj1 = MultivariatePolynomial(t);
            end
            M1 = obj1.matrixData;
            M2 = obj2.matrixData;
            numMon1 = obj1.numMonomials();
            numMon2 = obj2.numMonomials();
            numMonRes = numMon1 * numMon2;
            R = zeros(numMonRes, obj1.numIndeterminates()+1);
            for i=1:numMon1
                for j=1:numMon2
                    monIJ = obj1.matrixData(i, 1:end-1) + obj2.matrixData(j, 1:end-1);
                    coeffIJ = obj1.matrixData(i, end) * obj2.matrixData(j, end);
                    R((j-1)*numMon1 + i,:) = [monIJ coeffIJ];
                end
            end
            r = MultivariatePolynomial(R, obj1.ord, obj2.varNames);
        end % END mtimes
        
        function y = mpower(obj, b)            
            if b == 0                
                y = 1;
                return;
            end
            y = obj;
            for i = 1:b-1
                y =  obj * y;
            end
        end % END mpower
        
        
        function [q, r] = euclideanDivision(obj, divisors)
            numDiv = numel(divisors);
            r = MultivariatePolynomial(zeros(1, obj.numIndeterminates+1), obj.ord, obj.varNames);
            q = cell(numDiv, 1);
            for i=1:numDiv
                q{i} = MultivariatePolynomial(zeros(1, obj.numIndeterminates+1), obj.ord, obj.varNames);
            end
            p = obj;
            while ~p.iszzero()
                i = 1;
                divisionOccured = false;
                while i <= numDiv && ~divisionOccured
                    d = p.divide_lead_term_by_lead_term_of(divisors{i});
                    if ~isempty(d)
                        q{i} = q{i} + d;
                        p = p - d * divisors{i};
                        divisionOccured = true;
                    else
                        i = i + 1;
                    end
                end
                if ~divisionOccured
                    ltp = MultivariatePolynomial(p.leadTerm);
                    r = r + ltp;
                    p = p - ltp;
                end
            end
        end % END euclideanDivision
        
        function disp(obj)
            disp(['Polynomial of ', num2str(size(obj.matrixData,2)-1), ' variables'])
            if obj.iszzero
                disp('Zero polynomial')
            end
            polyString = '';
            for i=1:size(obj.matrixData, 1)
                coeff_i = obj.matrixData(i, end);
                monom_i = obj.matrixData(i, 1:end-1);
                sg = sign(coeff_i);
                monomStr_i = num2str(abs(coeff_i));
                if abs(coeff_i)==1 && sum(monom_i)>0
                    monomStr_i = '';
                end
                for j = 1:numel(monom_i)
                    if monom_i(j) >= 1
                        monomStr_i = [monomStr_i , obj.varNames{j}];
                        if monom_i(j) >= 2
                            monomStr_i = [monomStr_i, '^', num2str(monom_i(j))];
                        end
                    end
                    
                end
                if i==1 && sg==1
                    polyString = [polyString, monomStr_i, ''];
                else
                    if sg==1
                        polyString = [polyString, ' + ', monomStr_i, ''];
                    elseif sg==-1
                        polyString = [polyString, ' - ', monomStr_i, ''];
                    end
                end
            end
            disp(polyString )
        end
    end
end

