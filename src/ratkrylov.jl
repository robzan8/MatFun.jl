#=
canonical_cplx checks if the poles xi are ordered canonically.
=#
function canonical_cplx(xi::Vector{C})::Bool where {C<:Complex}
	m = length(xi)
	j = 1
	while j <= m
		if isreal(xi[j]) || isinf(xi[j])
			j = j+1
		else
			if j == m || (j < m && xi[j+1] != conj(xi[j]))
				return false
			end
			j = j+2
		end
	end
	return true
end

function [xi, mu, rho, eta] = poles_to_moebius(xi)
	% POLES_TO_MOEBIUS Moebius transformation with poles xi.
	%
	% Nonzero xi is replaced with (xi, mu) := (xi, 1) and
	% (rho, eta) := (1, 0), and xi = 0 is replaced by
	% (xi, mu) := (0, inf) and (rho, eta) := (0, 1).
	mu = ones(size(xi));
	rho = ones(size(xi));
	eta = zeros(size(xi));
	%rho = randn(size(xi));
	%eta = randn(size(xi));

	mu(xi == 0) = inf;
	xi(xi == 0) = 1;

	eta(isinf(mu)) = 1;
	rho(isinf(mu)) = 0;
end

% rat_krylov(A, b, xi, flag)
B = speye(N);
b = varargin{1}{2};
xi = varargin{1}{3};
m = size(xi, 2);
H = zeros(m+1, m);
K = zeros(m+1, m);
inner_product = @(x, y) y'*x;
if strcmp('real', varargin{1}{4}) realopt = 1;
else realopt = 0; end
rerun = 0;

U  = zeros(m+1, m);

% Cannot use real option if the poles are not ordered canonically.
if realopt && ~canonical_cplx(xi)
	realopt = 0;
	warning;
end

[xi, mu, rho, eta] = poles_to_moebius(xi);

V = zeros(N, m+1);

function run_krylov
	bd = false;
	bd_tol = eps(1);
	% Starting vector.
	V(:, 1) = b/induced_norm(b);
	j = 1;
	while j <= m
		% Computing the continuation combination.
		if j == 1, U(1, 1) = 1; else
			[Q, ~] = qr(K(1:j, 1:j-1)/mu(j) - H(1:j, 1:j-1)/xi(j));
			U(1:j, j) = Q(:, end);
		end

		% Compute new vector.
		w = V(:, 1:j)*U(1:j, j);
		w = rho(j)*(A*w) - eta(j)*(B*w);
		w = linear_solver(B/mu(j)-A/xi(j), w);

		% Orthogonalization.
		if isreal(xi(j)/mu(j)) || realopt == 0
			% MGS
			for reo = 0:1
				for reo_i = 1:j
					hh(1) = inner_product(w, V(:, reo_i));
					w = w - V(:, reo_i)*hh(1);
					H(reo_i, j) = H(reo_i, j) + hh(1);
				end
			end
			H(j+1, j) = induced_norm(w);
			V(:, j+1) = w/H(j+1, j);

			if abs(H(j+1, j)) < bd_tol*norm(H(1:j, j)), bd = true; break; end

			% Setting the decomposition.
			K(1:j+1, j) = H(1:j+1, j);
			K(1:j+1, j) = rho(j)*[U(1:j, j); 0] + K(1:j+1, j)/xi(j);
			H(1:j+1, j) = eta(j)*[U(1:j, j); 0] + H(1:j+1, j)/mu(j);
		else
			V(:, j+1) = real(w);
			V(:, j+2) = imag(w);

			% MGS
			for j = j:j+1
				for reo = 0:1
					for reo_i = 1:j
						hh(1) = inner_product(V(:, j+1), V(:, reo_i));
						V(:, j+1) = V(:, j+1) - V(:, reo_i)*hh(1);
						H(reo_i, j) = H(reo_i, j) + hh(1);
					end
				end
				H(j+1, j) = induced_norm(V(:, j+1));
				V(:, j+1) = V(:, j+1)/H(j+1, j);
				if abs(H(j+1, j)) < bd_tol*norm(H(1:j, j)), bd = true; break; end
			end

			% Setting the decomposition.
			rxi = real(1/xi(j-1)); ixi = imag(1/xi(j-1));
			cxi = [rxi ixi; -ixi rxi];

			rcnt = [real(U(1:j-1, j-1)) imag(U(1:j-1, j-1)); 0 0; 0 0];
			K(1:j+1, j-1:j) = H(1:j+1, j-1:j);
			K(1:j+1, j-1:j) = K(1:j+1, j-1:j)*cxi + rho(j-1)*rcnt;
			H(1:j+1, j-1:j) = H(1:j+1, j-1:j)/mu(j-1) + eta(j-1)*rcnt;
		end % realopt

		j = j+1;
	end % while j <= m

	if bd == true
		warning(['rat_krylov: ''lucky breakdown'' occured at iteration ' num2str(j)]);
		V = V(:, 1:j);
		K = K(1:j, 1:j-1);
		H = H(1:j, 1:j-1);
	end

end % run_krylov
