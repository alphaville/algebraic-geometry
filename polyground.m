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
assert(check.iszzero)
 

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
assert(check.iszzero)

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