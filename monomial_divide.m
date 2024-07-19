function y = monomial_divide(m1, m2)
y = m1 - m2;
if any(y<0)
    y = [];
end