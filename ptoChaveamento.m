function [y_1, v_1, ltx]= ptoChaveamento(vals)

syms v_1 y_1 v_0 y_0 v_2 y_2 a1 a2;
vars = [v_0, y_0, v_2, y_2, a1, a2];

[y_1, v_1] = solve(v_1^2 == v_0^2 + 2*a1*(y_1 - y_0), v_2^2 == v_1^2 + 2*a2*(y_2 - y_1), 'y_1','v_1');

ltx.y_1 = latex(simplify(y_1)); 
ltx.v_1 = latex(simplify(v_1));

y_1 = double(vpa(subs(y_1, vars, vals)));
v_1 = double(vpa(subs(v_1, vars, vals)));

end