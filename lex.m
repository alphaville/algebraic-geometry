function y = lex(a, b)
% Lexicographic ordering of multiindices
%
% INPUTS:
% a, b: multiindices (arrays)
%
% OUTPUT:
% 1 if a > b, -1 if b > a, 0 if a = b

d = a - b;
idx_d = find(d~=0, 1, 'first');
if isempty(idx_d)
    y = 0;
else
    y = sign(d(idx_d));
end