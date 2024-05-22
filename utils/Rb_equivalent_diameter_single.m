function Rb = Rb_equivalent_diameter_single(lambda_b, rb, rp, D)
% see eqn in Table 1 of Lai et al 2015
% original eqn comes from Gu et al 1998
% called equivalent diameter for single U-tube

% should use D or Ls in formula ???

Rb = 1 / (2*pi()*lambda_b) * log( rb/rp * sqrt(rp/D) );

end

