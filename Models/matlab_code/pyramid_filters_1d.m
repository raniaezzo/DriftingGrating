rng(0)

for input_type = {'white' 'pink' 'narrowband'}
    % How many time points? (power of 2)
    n = 10;
    T= 2^n;
    t = (0:T-1)-T/2;
    f = 0:T-1;

    % Which frequencies for the filters?
    f0 = 2.^[3 5 7];
    
    clear input_signal

    % Input signal contains unit sine waves at the frequencies of the filters
    switch input_type{1}
        case 'white'
            for ii = 1:length(f)/2
                input_signal(ii,:) = sin(2*pi*t/T * f(ii)+rand*2*pi);
            end
            input_signal = sum(input_signal)';
            
        case 'pink'
            for ii = 1:length(f)/2
                input_signal(ii,:) = 1/f(ii+1)*sin(2*pi*t/T * f(ii+1)+rand*2*pi);
            end
            input_signal = sum(input_signal)';
        case 'narrowband'
            input_signal = sum(sin(2*pi*t/T .* f0'))';
    end

    % Generate Gabor filters in quadrature pairs
    for filter_normalization = {'Equal Height' 'Equal Area'}

        for ii = 1:length(f0)

            sz = T/2/f0(ii);
            sine_wave         = sin(2*pi*t/length(t)* f0(ii));
            cosine_wave       = cos(2*pi*t/length(t)* f0(ii));
            gaussian_enevelope  = exp(-1/2*(t/sz).^2);

            switch filter_normalization{1}
                case 'Equal Height'
                    % do nothing
                case 'Equal Area'
                    gaussian_enevelope = gaussian_enevelope / sum(gaussian_enevelope);
            end

            % gaussian_enevelope to make it spatially specific
            sine_filter(:,ii)   = gaussian_enevelope.*sine_wave;
            cosine_filter(:,ii) = gaussian_enevelope.*cosine_wave;

        end

        for ii = 1:length(f0)
            output_signals(ii,:, :) = ...
                [conv(input_signal, sine_filter(:,ii), 'same') ...
                conv(input_signal, cosine_filter(:,ii), 'same')];
        end
        output_energy = sum(output_signals.^2,3)';

        % Plotting parameters
        lw = 1.5;  % linewidth
        fs = 12; % fontsize

        figure, h = tiledlayout(3,2);
        title(h, sprintf('Filters: %s\n Input: %s', ...
            filter_normalization{1}, input_type{1}));

        % Input Signal (time domain)
        nexttile();
        plot(t/T, input_signal, LineWidth=lw, Color='k')
        title('Input Signal (time domain)')
        xlabel('Time (s)')
        ylabel('Intensity')
        set(gca, 'FontSize', fs)

        % Input Signal (frequency domain)
        nexttile();
        plot(f, abs(fft(input_signal)), LineWidth=lw, Color='k')
        title('Input Signal (freq domain)')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude')
        set(gca, 'FontSize', fs)
        xlim([0 max(f)/2])

        % Filter (time domain)
        nexttile();
        plot(t/T, sine_filter+ 2*max(sine_filter(:))*(1:length(f0)), LineWidth=lw)
        %plot(t/T, sine_filter+ 2*max(sine_filter(:))*(1:length(f0)), LineWidth=lw)
        %hold on
        %plot(t/T, cosine_filter+ 2*max(sine_filter(:))*(1:length(f0)), LineWidth=lw)
        xlim([-0.25, 0.25])
        title('Filters (time domain)')
        xlabel('Time (s)')
        ylabel('Response')
        set(gca, 'FontSize', fs)

        % Filter (frequency domain)
        nexttile();
        f = 0:T-1;
        plot(f, abs(fft(sine_filter)), LineWidth=lw);
        xlim([0 max(f)/2])
        title('Filters (freq domain)')
        xlabel('Frequency (Hz)')
        ylabel('Response')
        set(gca, 'FontSize', fs)


        % Output (time domain)
        nexttile();
        plot(t, output_energy, LineWidth=lw)
        title('Output (time domain)')
        xlabel('Time (s)')
        ylabel('Energy')
        set(gca, 'FontSize', fs)

        % Output (frequency domain)
        nexttile();
        plot(f, abs(fft(output_signals(:,:,1), [], 2)), LineWidth=lw)
        title('Output (freq domain)')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude')
        set(gca, 'FontSize', fs)
        xlim([0 max(f)/2])


    end
end