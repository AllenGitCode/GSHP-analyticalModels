function G = Gfunction_FLS(lambda_g, alpha_g, H, rb, z, t)
% can use z = H/2 to get the average G-function value

fun = @(zprime, rb, z, t, alpha_g) 1./(sqrt( rb.^2 + (z-zprime).^2 )) .* erfc( (sqrt( rb.^2 + (z-zprime).^2 ))./sqrt(4*alpha_g*t) );

q1 = integral(@(zprime) fun(zprime, rb, z, t, alpha_g), 0, H, 'ArrayValued',true);
q2 = integral(@(zprime) fun(zprime, rb, z, t, alpha_g), -H, 0, 'ArrayValued',true);

G = 1 / (4*pi()*lambda_g) * (q1 - q2);

% examples:
%   Gfunction_FLS(lambda_g, alpha_g, H, 0.5, [1 2 3 4 5 6 7 8 9 10]', 24*60*60)
%   Gfunction_FLS(lambda_g, alpha_g, H, 0.5, H/2, 24*60*60)
end

