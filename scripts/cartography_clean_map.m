function cartography_clean_map(input_file, output_file, cartography_root)
% I/O and coefficient permutation only. The sky synthesis is upstream makemap.
% Adapted from maps_conversion_nocleandistribution.m's clean-map branch.
matlab_root = fullfile(cartography_root, 'complex_analysis', 'matlab');
addpath(matlab_root);
addpath(fullfile(matlab_root, 'Utilities'));
addpath(fullfile(matlab_root, 'spherelib'));
input = load(input_file);
clms = input.clms(:);
lmax = sqrt(length(clms)) - 1;
assert(lmax == floor(lmax));
[lvec, mvec] = getLMvec(lmax);
lvec_pta = []; mvec_pta = [];
for ll=0:lmax
    for mm=-ll:ll
        lvec_pta = [lvec_pta ll];
        mvec_pta = [mvec_pta mm];
    end
end
pOpt = 0*clms;
permutation = zeros(length(clms),1);
for ii=1:length(clms)
    idx = find(lvec_pta==lvec(ii) & mvec_pta==mvec(ii));
    permutation(ii) = idx;
    pOpt(ii) = clms(idx);
end
% plm2xyz writes a recalculable Legendre cache relative to the working directory.
if ~exist('LEGENDRE', 'dir'), mkdir('LEGENDRE'); end
% Preserve the monopole and signed values. No uncertainty or S/N calculations.
[clean_normalized, RA_hours, DEC_deg, dOmega] = makemap(pOpt, 1);
save('-v7', output_file, 'clean_normalized', 'RA_hours', 'DEC_deg', ...
     'dOmega', 'pOpt', 'permutation', 'lvec', 'mvec');
end
