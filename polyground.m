M = [4     1     0     1
     3     3     5     2
     3     3     5    -4
     3     3     5     8
     2     1     1     5
     2     0     1     4
     2     0     1     3
     2     0     1    -1
     1     1     1     7
     1     1     1    10
     1     1     1    1
     1     0     0    5
     0     0     0    4];
 
ord = @(s1, s2) grlex(s1, s2);
p = MultivariatePolynomial(M, ord);

M2 = [1     0     0    56
      0     0     1    4
      0     0     0   -3];
q = MultivariatePolynomial(M2, ord);

% s = p + q;
% s.matrixData

w = p * q;

%% Exercise 1
clc; clear;
ord = @(s1, s2) grlex(s1, s2);
f = MultivariatePolynomial([7 2 1;  3 2 1; 0 1 -1; 0 0 1], ord, {'x', 'y'});
f1 = MultivariatePolynomial([1 2 1; 1 0 -1], ord, {'x', 'y'});
f2 = MultivariatePolynomial([1 0 1; 0 3 -1], ord, {'x', 'y'});

divisors = {f2, f1};
[q, r] = f.euclideanDivision(divisors);
q{1}
q{2}
r

 
check = q{1}*divisors{1} + q{2}*divisors{2} + r - f;
assert(check.iszero)
 

%% Exercise 2
clc; clear;
ord = @(s1, s2) lex(s1, s2);
f = MultivariatePolynomial([1 2 2 1; 1 1 0 1; 0 1 1 -1], ord, {'x', 'y', 'z'});
f1 = MultivariatePolynomial([1 0 0 1; 0 2 0 -1], ord, {'x', 'y', 'z'});
f2 = MultivariatePolynomial([0 1 0 1; 0 0 3 -1], ord, {'x', 'y', 'z'});
f3 = MultivariatePolynomial([0 0 2 1; 0 0 0 -1], ord, {'x', 'y', 'z'});

divisors = {f2, f3, f1};
[q, r] = f.euclideanDivision(divisors);


q{1}
q{2}
q{3}
r

check = q{1}*divisors{1} + q{2}*divisors{2} + q{3}*divisors{3} + r - f;
assert(check.iszero)

%% Demo
clc; clear;
ord = @(s1, s2) grlex(s1, s2);
f = MultivariatePolynomial([2 5 2; 4 3 -1; 5 6 8], ord, {'x', 'y'});
f1 = MultivariatePolynomial([2 4 1], ord, {'x', 'y'});
f2 = MultivariatePolynomial([3 2 1], ord, {'x', 'y'});
f3 = MultivariatePolynomial([5 1 1], ord, {'x', 'y'});


divisors = {f1, f2, f3};
[q, r] = f.euclideanDivision(divisors);
q{:}
r

%% Exercise 1, sec 6, p 81
clc
ord = @(s1, s2) lex(s1, s2);

[x, y, z] = polysymbols({'x', 'y', 'z'}, ord);

g1 = x*y^2 - x*z + y;
g2 = x*y - z^2;
g3 = x - y*z^4;


g = g1 - y*g2 + z*g3;
lt_g = g.leadTermAsPolynomial
lt_gi  = {g1.leadTermAsPolynomial, g2.leadTermAsPolynomial, g3.leadTermAsPolynomial};

[q, r] = lt_g.euclideanDivision(lt_gi) ;
r

%%
clc
clear

ord = @(s1, s2) grlex(s1, s2);

vars = {'x', 'y', 'z'};
[x, y, z] = polysymbols(vars, ord);
id = MultivariatePolynomial.id(vars, ord);

f1 = 3*x - z^4;
f2 = y - z^5;
f3 = y - z^5;

F = {f1, f2, f3};
I = PolynomialIdeal(F);
I.grobnerBasis(true);
G = I.grobner;
G{:}

%I.ismember(x^2 - x^5*z^5 - x*y^3*z^4 + x^5*y + x^2*y^3 - x*z^4 + 2*z^4 - 2*x)



%%
clc
clear

ord = @(s1, s2) lex(s1, s2);
vars = {'x', 'y', 'z'};
[x, y, z] = polysymbols(vars, ord);
id = MultivariatePolynomial.id(vars, ord);

f1 = x^6 + 3*x^4*y^2 + 3*x^2*y^4 - 4*x^2*y^2 + y^6;
f2 = 16*x^4 - 8*x^2 - x*y^2 - 5*y^3 + 1.5*y^2 + id;

