function varargout = polysymbols(varNames, ord)

if nargin==1
    ord = @(s1, s2) lex(s1, s2);
end

n = numel(varNames);
varargout = cell(n, 1);
for i = 1:n
    d = zeros(1, n+1);
    d(i) = 1; d(end) = 1;
    varargout{i} = MultivariatePolynomial(d, ord, varNames);
end
