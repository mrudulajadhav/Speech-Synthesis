    % read the input audio
    [A, Fs]= audioread('phrase.wav');
    A = A(:,1);
    frame_size = 0.030; % 30 ms
    frame_length = Fs*frame_size;

    percent_overlap = 0; % 50% overlap
    frame_overlap = frame_length*percent_overlap/100;
    frame_step = frame_length - frame_overlap;
    audio_length = length(A);

    no_of_frames = floor(abs(audio_length - frame_overlap)/(abs(frame_length - frame_overlap)));
    rest_samples = mod(abs(audio_length - frame_overlap) , abs(frame_length - frame_overlap));

    % padding if remaining samples != 0
    if rest_samples ~= 0
             pad_signal_length = int16(frame_length - rest_samples);
             z = zeros(pad_signal_length,1);
             %A(end+1:end+pad_signal_length)= 0;
             A = vertcat(A, zeros(pad_signal_length,1));
             no_of_frames = no_of_frames + 1;
    end

    % framing
    mat = zeros(no_of_frames, frame_length);
    i = 1;

    for r = 1:no_of_frames
        for c = 1:frame_length
            mat(r,c) = A(i);
            i = i+1;
        end
        i = i-frame_overlap;
    end

    % hamming window coefficients
    hamming = zeros(1,frame_length);
    for i=1:frame_length
        x = 2*pi*(i-1);
        x = x./(frame_length-1);
        hamming(1,i)=0.54-0.46.*cos(x);
    end

    % windowing
    for r = 1:no_of_frames
        mat(r,:)= mat(r,:).*hamming;
    end
 
    % computing ste and ZCR
    ZCR = zeros(no_of_frames,1);
    STE = zeros(no_of_frames,1);
    for r = 1:no_of_frames
        STE(r) = sum(mat(r,:).^2);
        ZCR(r) = mean(abs(diff(sign(mat(r,:)))));
    end
    
    STE_thresh = 0.456;
    ZCR_thresh = 0.5;

    OUTPUT_STE = zeros(no_of_frames,1);
    for i = 1:no_of_frames
        if STE(i) > STE_thresh
            OUTPUT_STE(i) = 1.0;
        else
            OUTPUT_STE(i) = 0.0;
        end
    end
    
    OUTPUT_ZCR = zeros(no_of_frames,1);
    for i = 1:no_of_frames
        if ZCR(i) > ZCR_thresh
            OUTPUT_ZCR(i) = 0.0;
        else
            OUTPUT_ZCR(i) = 1.0;
        end
    end
   
    % output array for plot
    o_arr1 = zeros(length(A),1);
    j = 1;k = 1;
    for i = 1:no_of_frames
        while k <= frame_length
            if OUTPUT_STE(i) == 1.0
                o_arr1(j) = 0.2;
            else
                o_arr1(j) = -0.2;
            end
            k = k + 1;
            j = j + 1;
        end
        k = 1 + frame_overlap;
    end
    
    % output array for plot
    o_arr2 = zeros(length(A),1);
    j = 1;k = 1;
    for i = 1:no_of_frames
        while k <= frame_length
            if OUTPUT_ZCR(i) == 1.0
                o_arr2(j) = 0.3;
            else
                o_arr2(j) = -0.3;
            end
            k = k + 1;
            j = j + 1;
        end
        k = 1 + frame_overlap;
    end
 
    x = linspace(1,length(A),length(A));
    figure(1);
    plot(x,A)
    hold on
    plot(x,o_arr1,'-r')
    plot(x,o_arr2,'-g')
    legend('INPUT SIGNAL','OUTPUT WITH STE','OUTPUT WITH ZCR');
    title('VOICED UNVOICED SEGMENTS')
    hold off
   
    
    % lpc for each frame (row wise)
    no_of_coeff = 13; %no of lpc coeff
    lpc_mat = zeros(no_of_frames, no_of_coeff);
    for i = 1:no_of_frames
       [lpc_coeff,g] = lpc(mat(i,:),no_of_coeff); 
       lpc_mat(i,:) = lpc_coeff(2:end);
    end
    
    % excitation for voiced
    %creating excitation (impulse train with 7.5ms pitch)
    exciteV = zeros(1,frame_length); % (framelength)
    fs = Fs; %sampling freq is 8kHz
    pitch_period = 0.0075; %7.5 ms
    no_of_samples_in_one_period = fs*pitch_period;
    oneperiod = randi([0 1], 1,no_of_samples_in_one_period); %generating impulse of samples
    p = 0;

    l = frame_length/no_of_samples_in_one_period;
    for i = 1:l
        for j = 1:no_of_samples_in_one_period
            exciteV(1,p+j) = oneperiod(1,j);
        end
        p = p+no_of_samples_in_one_period;
    end

    % excitation for unvoiced
    exciteUV = normrnd(0,1,[1,frame_length]);

    synth_signal = zeros(1,length(A));
    p=1;
    for i = 1:no_of_frames
        if OUTPUT_ZCR(i)==1
           s = filter(1,[1 lpc_mat(i,:)], exciteV);
        else
           s = filter(1,[1 lpc_mat(i,:)], exciteUV);
        end
        synth_signal(p:p+frame_length-1)= s;
        p = p+frame_length;
    end

    %plot synthesised and original signal
    x = linspace(1,length(A),length(A));
    
    figure(2);
    subplot(2,1,1);
    plot(x,A);
    title('SIGNAL PLOT')
    legend('original signal');
    subplot(2,1,2);
    plot(x,synth_signal,'-r');
    legend('synthesized signal');

    %plot difference signal
    A = reshape(A,[1,length(A)]);  
    rmse = sqrt(mean((A - synth_signal).^2));
    figure(3);
    plot(x, A-synth_signal);
    title('DIFFERENCE SIGNAL')
    
    synth_signal = rescale(synth_signal,-1,1);    
    audiowrite('synth_signal_low_pitch.wav',synth_signal,8000);
    soundsc(synth_signal);
