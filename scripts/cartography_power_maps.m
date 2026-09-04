function cartography_power_maps(input_file, output_file, cartography_root)
% Convert PTA-order clean coefficients, dirty-map vector, and Fisher matrix
% to the upstream 1-degree clean map and one requested radiometer scalar.
% This function never constructs or saves a radiometer sky map.
% Adapted from maps_conversion_nocleandistribution.m without S/N or trials.
matlab_root = fullfile(cartography_root, 'complex_analysis', 'matlab');
addpath(matlab_root);
addpath(fullfile(matlab_root, 'Utilities'));
addpath(fullfile(matlab_root, 'spherelib'));
input = load(input_file);
clms_pta = input.clms(:);
X_pta = input.X(:);
M_pta = input.M;
target_RA_hours = input.target_RA_hours;
target_DEC_deg = input.target_DEC_deg;
lmax = sqrt(length(clms_pta)) - 1;
assert(lmax == floor(lmax));
assert(all(size(M_pta) == [length(clms_pta), length(clms_pta)]));
assert(length(X_pta) == length(clms_pta));

[lvec, mvec] = getLMvec(lmax);
lvec_pta = [];
mvec_pta = [];
for ll=0:lmax
    for mm=-ll:ll
        lvec_pta = [lvec_pta ll];
        mvec_pta = [mvec_pta mm];
    end
end

% Reorder from PTA l-major, m=-l,...,+l to the LIGO convention used by
% makemap and diagPixel.
permutation = zeros(length(clms_pta), 1);
for ii=1:length(clms_pta)
    permutation(ii) = find(lvec_pta==lvec(ii) & mvec_pta==mvec(ii));
end
clms = clms_pta(permutation);
X = X_pta(permutation);
M = M_pta(permutation, permutation);

% plm2xyz writes a recalculable Legendre cache relative to the working
% directory.
legendre_cache = fullfile(pwd, 'LEGENDRE');
if ~exist(legendre_cache, 'dir'), mkdir(legendre_cache); end
[clean_normalized, RA_hours, DEC_deg, dOmega] = makemap(clms, 1);
target_pixel = find(abs(RA_hours-target_RA_hours) < 1e-12 ...
                  & abs(DEC_deg-target_DEC_deg) < 1e-12);
assert(length(target_pixel) == 1);

% Construct only the requested row of the pixel conversion matrix.  A full
% 65160x49 complex matrix exceeds this pinned Octave build's reliable memory
% path, while the radiometer value requested by the notebook needs one row.
target_U = zeros(1, length(clms));
for ii=1:length(clms)
    basis = zeros(length(clms), 1);
    basis(ii) = 1;
    basis_map = makemap(basis, 1, 0, 1);
    target_U(ii) = basis_map(target_pixel);
end

% This is diagPixel's U*M*U' definition written elementwise.  Complex
% row-by-column BLAS products are broken in the pinned Octave 9.4 build,
% while the equivalent row-times-matrix and elementwise sum are stable.
target_UM = target_U*M;
target_diag_fisher = real(sum(target_UM.*conj(target_U)));
target_X_pixel = sum(target_U.*transpose(X));
assert(target_diag_fisher > 0);
target_radiometer_normalized = real(target_X_pixel./target_diag_fisher);
target_radiometer_sigma_normalized = real(target_diag_fisher.^-0.5);

save('-v7', output_file, 'clean_normalized', 'target_radiometer_normalized', ...
     'target_radiometer_sigma_normalized', 'target_diag_fisher', ...
     'target_U', 'target_RA_hours', 'target_DEC_deg', ...
     'target_pixel', 'RA_hours', ...
     'DEC_deg', 'dOmega', 'clms', 'X', 'M', 'permutation', 'lvec', 'mvec');
end
