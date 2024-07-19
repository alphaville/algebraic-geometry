function x = monomial_sort(x, ord)
% INPUTS:
% x - matrix of multiindices (monomials) and coefficients, where the LAST
%     column stores the coefficients
% ord - monomial ordering as function handle
%       defaults to lex
% OUTPUT:
% x - sorted x

if nargin == 1
    ord = @(s1, s2) lex(s1, s2);
end

% `beginIdx` and `endIdx` marks the first and last index to check
beginIdx = 1;
endIdx = size(x, 1) - 1;

while beginIdx <= endIdx
    newBeginIdx = endIdx;
    newEndIdx = beginIdx;
    for ii = beginIdx:endIdx
        if ord(x(ii,1:end-1), x(ii + 1, 1:end-1)) == -1
            [x(ii+1, :), x(ii, :)] = deal(x(ii, :), x(ii+1, :));
            newEndIdx = ii;
        end
    end
    
    % decreases `endIdx` because the elements after `newEndIdx` are in correct order
    endIdx = newEndIdx - 1;
    
    for ii = endIdx:-1:beginIdx
        if ord(x(ii, 1:end-1), x(ii + 1,1:end-1)) == -1
            [x(ii+1, :), x(ii, :)] = deal(x(ii, :), x(ii+1,:));
            newBeginIdx = ii;
        end
    end
    % increases `beginIdx` because the elements before `newBeginIdx` are in correct order
    beginIdx = newBeginIdx + 1;
end
end