function norms = colNorm(H)

norms = zeros(cols(H), 1);

for i = 1:cols(H)
  norms(i) = norm(H(:, i));
end
