% Parameters
N = 768;
numscales = 6;
numorientations = 4;

[xx, yy] = meshgrid(0:N-1, 0:N-1);

% spatial frequencies we will probe (cycles across full image width)
test_cycles = [12 24 48 96 192];  % you can extend this
nTests = numel(test_cycles);

band_energy_all = zeros(numscales, nTests);

for tc = 1:nTests

    nCyclesAcross = test_cycles(tc);

    % convert that to cycles/pixel magnitude
    f_cycPerPix = nCyclesAcross / N;

    % make a horizontal grating (orientation doesn't actually matter for radial freq)
    % horizontal bars = variation along y, so kx=0, ky=f
    kx = 0;
    ky = f_cycPerPix;
    stim = cos(2*pi*(kx*xx + ky*yy));  % N x N

    % build pyramid
    coeff_all = buildSCFpyr(stim, numscales+2, numorientations);

    % compute total energy per band, summed across orientations
    for ii = 2:(numscales+1)          % ii = 2..7
        band_s = coeff_all{ii};       % 1 x numorientations
        subcat = cat(3, band_s{:});   % stack orientations
        e = sum(abs(subcat(:)).^2);   % squared magnitude energy
        band_energy_all(ii-1, tc) = e;
    end
end

% normalize each column so we can compare peaks
band_energy_norm = band_energy_all ./ max(band_energy_all, [], 1);

% print which band wins for each tested frequency
fprintf('Tested frequency (cycles/image) -> best band (ii index in coeff_all):\n');
for tc = 1:nTests
    [~, best_band] = max(band_energy_all(:, tc));
    fprintf('%5d cycles -> band %d (coeff_all{%d})\n', ...
        test_cycles(tc), best_band, best_band+1);
end

% % optional: visualize tuning curves
% figure; hold on;
% colors = lines(numscales);
% for b = 1:numscales
%     plot(test_cycles, band_energy_norm(b,:), '-o', 'LineWidth', 2, 'Color', colors(b,:));
% end
% xlabel('Input grating frequency (cycles across full image)');
% ylabel('Normalized band energy');
% title('Empirical spatial frequency tuning per pyramid band');
% legend(arrayfun(@(b) sprintf('band %d (coeff\\\\_all{%d})', b, b+1), 1:numscales, 'UniformOutput', false), ...
%        'Location', 'NorthWest');
% grid on;
