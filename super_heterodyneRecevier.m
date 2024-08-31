[x, fs_x] = audioread('Short_RussianVoice.wav');
[y, fs_y] = audioread('Short_SkyNewsArabia.wav');

%--------------------trasmitter(AM modulation)---------------------------%
% Set carrier frequencies for FDM
fc_x = 100e3; % 100 kHz
fc_y = 155e3; % 155 kHz

% Increase sampling rate 
upsample_factor = 16;

% Determine the shorter and longer signals
if length(x) <= length(y) % skynews is longer
    shorter_signal = x;
    longer_signal = y;
    fs_shorter = fs_x;
    fs_longer = fs_y;
    % Pad the shorter signal with zeros to match the length of the longer one
    shorter_signal = [shorter_signal; zeros(length(longer_signal) - length(shorter_signal), size(shorter_signal, 2))];
   % Check if carrier frequency for x exceeds sampling limit
   if fc_x > fs_x / 2 || fc_y > fs_y / 2
    % channel stream
    shorter_signal = shorter_signal(:,1) + shorter_signal(:,2);
    longer_signal = longer_signal(:,1) + longer_signal(:,2);
   
    % Interpolate signal x
    shorter_signal = interp(shorter_signal, upsample_factor);
    longer_signal = interp(longer_signal, upsample_factor);
    
    % Update the length and sampling rate
    fs_shorter = fs_shorter * upsample_factor;
    fs_longer = fs_longer * upsample_factor;
    
    % Generate carriers for FDM
    t_x = (0:length(longer_signal)-1)/fs_longer; % Use the longer signal's sampling rate
    t_y = (0:length(longer_signal)-1)/fs_longer; % Use the longer signal's sampling rate
    carrier_x = cos(2*pi*fc_x*t_x);
    carrier_y = cos(2*pi*fc_y*t_y);
    % Amplitude Modulation for FDM
    modulated_x = shorter_signal .* carrier_x(:);
    modulated_y = longer_signal.* carrier_y(:);
   end

   else % russian is longer
    shorter_signal = y;
    longer_signal = x;
    fs_shorter = fs_y;
    fs_longer = fs_x;
   % Pad the shorter signal with zeros to match the length of the longer one
    shorter_signal = [shorter_signal; zeros(length(longer_signal) - length(shorter_signal), size(shorter_signal, 2))];

    % Check if carrier frequency for y or x exceeds sampling rate
    if fc_x > fs_x / 2 || fc_y > fs_y / 2 

    % channel stream
    shorter_signal = shorter_signal(:,1) + shorter_signal(:,2);
    longer_signal = longer_signal(:,1) + longer_signal(:,2);
    
    % Interpolate signal y
    longer_signal = interp(longer_signal, upsample_factor_y);
    shorter_signal = interp(shorter_signal, upsample_factor_x);
    
    % Update the length and sampling rate
    fs_longer = fs_longer * upsample_factor_y;
    fs_shorter = fs_shorter * upsample_factor_x;

    % Generate carriers for FDM
    t_x = (0:length(longer_signal)-1)/fs_longer; % Use the longer signal's sampling rate
    t_y = (0:length(longer_signal)-1)/fs_longer; % Use the longer signal's sampling rate
    carrier_x = cos(2*pi*fc_x*t_x);
    carrier_y = cos(2*pi*fc_y*t_y);
    % Amplitude Modulation for FDM
    modulated_x = longer_signal .* carrier_x(:);
    modulated_y = shorter_signal.* carrier_y(:);
   end
end

% Frequency Division Multiplexing (FDM)
fdm_signal = modulated_x + modulated_y;

% Plot the frequency domain representation of the FDM signal
figure;
% Define frequencies_fdm based on fs_longer
frequencies_fdm_range = linspace(-fs_longer/2, fs_longer/2, length(fdm_signal));
spectrum_fdm = fftshift(fft(fdm_signal));
plot(frequencies_fdm_range, spectrum_fdm);
title('FDM Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-350e3, 350e3]); % Set X-axis limits
ylim([0, 55e3]); % Set Y-axis limits
grid on;

%-----------------------------Design_For_BPF_Filters-------------------%
  % Filters design
filter_order = 100; 
center_frequency_x = 100e3; %  carrier at 100 kHz
center_frequency_y = 155e3; % carrier at 155 kHz
center_frequency_IF = 27.5e3; % carrier at 27.5 kHz
bandwidth = 20e3; % Bandwidth of the filter

% Design bandpass filters
filter_for_x = fir1(filter_order, [(center_frequency_x - bandwidth) (center_frequency_x + bandwidth)] / (fs_longer/2), 'bandpass');
filter_for_y = fir1(filter_order, [(center_frequency_y - bandwidth) (center_frequency_y + bandwidth)] / (fs_longer/2), 'bandpass');
filter_for_IF = fir1(filter_order, [(center_frequency_IF - bandwidth) (center_frequency_IF + bandwidth)] / (fs_longer/2), 'bandpass');

%------------------------------RECEIVER----------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
russian_channel = 1;
skynews_channel = 2;

% make the user choose between russian channel or SkyNewsArabia channel %
selected_channel = skynews_channel ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if selected_channel == russian_channel
      %-------------------------------RF_Stage---------------------------------%
      % % Apply BPF filter
      filtered_x = filter(filter_for_x, 1, fdm_signal);

      % Plot filtered signal
      figure;
      plot(frequencies_fdm_range, abs(fftshift(fft(filtered_x))));
      title('RUSSIAN at carrier frequency');
      xlabel('Frequency (Hz)');
      ylabel('Magnitude');
      xlim([-350e3, 350e3]); % Set X-axis limits
      ylim([0, 55e3]); % Set Y-axis limits
      grid on;
      %-------------------------------IF_Stage---------------------------------%
      % define the mixer
      if_freq= 27.5e3; % Adjusted IF frequency

      % Define the mixer signal
      mixer_x = cos(2*pi*(if_freq +fc_x)*t_x);

      % Down conversion to IF frequency
      down_x_IF = fdm_signal .* mixer_x';

      % Plot the converted signals at IF
      figure;
      frequencies_filtered_range_x = linspace(-fs_longer/2 , fs_longer/2, length(down_x_IF));
      plot(frequencies_filtered_range_x, abs(fftshift(fft(down_x_IF))));
      title('RUSSIAN at IF before BPF (NO RF)');
      xlabel('Frequency (Hz)');
      ylabel('Magnitude');
      xlim([-350e3, 350e3]); % Set X-axis limits
      ylim([0, 35e3]); % Set Y-axis limits
      grid on;

      % apply IF bandpass filter to reject the high frequency components
      % resulted from the down conversion operation
      filtered_Russian = filter(filter_for_IF, 1, down_x_IF);

      % Plot signals after down conversion to IF without the higher freq components
      figure;
      plot(frequencies_fdm_range, abs(fftshift(fft(filtered_Russian))));
      title('RUSSIAN at IF after BPF (NO RF) ');
      xlabel('Frequency (Hz)');
      ylabel('Magnitude');
      xlim([-350e3, 350e3]); % Set X-axis limits
      ylim([0, 33e3]); % Set Y-axis limits
      grid on;

      %-----------------------------BaseBand_Detection_Stage---------------------%
      %define the baseband mixer
      baseband_mixer =  cos(2*pi*(if_freq)*t_x);

      % down conversion to the baseband
      baseband_signal_x = filtered_Russian.*baseband_mixer';

      %plot signal at baseband
      figure;
      plot(frequencies_fdm_range, abs(fftshift(fft(baseband_signal_x))));
      title('RUSSIAN at baseband before the LPF (NO RF) ');
      xlabel('Frequency (Hz)');
      ylabel('Magnitude');
      xlim([-350e3, 350e3]); % Set X-axis limits
      ylim([0, 45e3]); % Set Y-axis limits
      grid on;
      % lowpass Filter design
      filter_order = 100; % Filter order
      cutoff_frequency = 25e3; % Cutoff frequency 
      LPfilter = fir1(filter_order, cutoff_frequency / (fs_longer/2), 'low');

      % apply the lowpass filter
      filtered_baseband_signal_x = filter(LPfilter , 1,baseband_signal_x );

      %plot signal at baseband
      figure;
      plot(frequencies_fdm_range, abs(fftshift(fft(filtered_baseband_signal_x))));
      title('RUSSIAN at baseband (NO RF)');
      xlabel('Frequency (Hz)');
      ylabel('Magnitude');
      xlim([-350e3, 350e3]); % Set X-axis limits
      ylim([0, 45e3]); % Set Y-axis limits
      grid on;


      %--------------------------hear sky news arabia channel--------------------% 
      filtered_baseband_signal_x = 8 .* filtered_baseband_signal_x; % gain
      filtered_baseband_signal_x = decimate(filtered_baseband_signal_x,16); % resample to origin
      fs_longer = fs_longer/16; % resample to origin
      sound(filtered_baseband_signal_x , fs_longer);


else % the chosen channel is skynews arabia %

      %-------------------------------RF_Stage---------------------------------%
          % Apply BPF filter
          filtered_y = filter(filter_for_y, 1, fdm_signal);

          % Plot filtered signal
          figure ;
          plot(frequencies_fdm_range, abs(fftshift(fft(filtered_y))));
          title('SKY NEWS at carrier frequency');
          xlabel('Frequency (Hz)');
          ylabel('Magnitude');
          xlim([-350e3, 350e3]); % Set X-axis limits
          ylim([0, 25e3]); % Set Y-axis limits
          grid on;

      %-------------------------------IF_Stage---------------------------------%
          % define the mixer
          if_freq= 27.5e3; % Adjusted IF frequency

          % Define the mixer signal
          mixer_y = cos(2*pi*(if_freq + fc_y)*t_y);

          % Down conversion to IF frequency
          down_y_IF = fdm_signal.* mixer_y';

          % Plot the converted signal at IF
          figure ;
          frequencies_filtered_range_y = linspace(-fs_longer/2, fs_longer/2, length(down_y_IF));
          plot(frequencies_filtered_range_y, abs(fftshift(fft(down_y_IF))));
          title('SKY NEWS ARABIA at IF before BPF (NO RF)');
          xlabel('Frequency (Hz)');
          ylabel('Magnitude');
          xlim([-400e3, 400e3]); % Set X-axis limits
          ylim([0, 15e3]); % Set Y-axis limits
          grid on;

          % apply IF bandpass filter to reject the high frequency components
          % resulted from the down conversion operation
          filtered_skyNews= filter(filter_for_IF, 1, down_y_IF);

          % Plot signal after down conversion to IF without the higher freq components
          figure ;
          plot(frequencies_fdm_range, abs(fftshift(fft(filtered_skyNews))));
          title('SKY NEWS ARABIA at IF after BPF (NO RF)');
          xlabel('Frequency (Hz)');
          ylabel('Magnitude');
          xlim([-400e3, 400e3]); % Set X-axis limits
          ylim([0, 15e3]); % Set Y-axis limits
          grid on;

      %-----------------------------BaseBand_Detection_Stage---------------------%
          % define the baseband mixer
          baseband_mixer =  cos(2*pi*(if_freq)*t_x);

          % down conversion to the baseband
          baseband_signal_y = filtered_skyNews.*baseband_mixer';

          %plot signal at baseband
          figure ;
          plot(frequencies_fdm_range, abs(fftshift(fft(baseband_signal_y))));
          title('SKY NEWS ARABIA at baseband before LPF (NO RF)');
          xlabel('Frequency (Hz)');
          ylabel('Magnitude');
          grid on;
          xlim([-400e3, 400e3]); % Set X-axis limits
          ylim([0, 15e3]); % Set Y-axis limits
          % lowpass Filter design
          filter_order = 100; % Filter order
          cutoff_frequency = 25e3; % Cutoff frequency 
          LPfilter = fir1(filter_order, cutoff_frequency / (fs_longer/2), 'low');

          %apply the lowpass filter
          filtered_baseband_signal_y = filter(LPfilter , 1,baseband_signal_y);

          %plot signal at baseband
          figure ;
          plot(frequencies_fdm_range, abs(fftshift(fft(filtered_baseband_signal_y))));
          title('SKY NEWS ARABIA at baseband (NO RF)');
          xlabel('Frequency (Hz)');
          ylabel('Magnitude');
          grid on;
          xlim([-400e3, 400e3]); % Set X-axis limits
          ylim([0, 15e3]); % Set Y-axis limits

      %--------------------------hear sky news arabia channel--------------------%
          filtered_baseband_signal_y = 8 .* filtered_baseband_signal_y; % gain
          filtered_baseband_signal_y = decimate(filtered_baseband_signal_y,16); % resample to origin
          fs_longer = fs_longer/16; % resample to origin
          sound(filtered_baseband_signal_y , fs_longer);

end % end of if condition statement %