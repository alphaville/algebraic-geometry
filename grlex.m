function y = grlex(a, b)
% Graded lexicographic ordering of multiindices
%
% INPUTS:
% a, b: multiindices (arrays)
%
% OUTPUT:
% 1 if a > b, -1 if b > a, 0 if a = b

d_deg = sum(a) - sum(b);
if d_deg ~= 0
    y = sign(d_deg);
else
    y = lex(a, b);
end