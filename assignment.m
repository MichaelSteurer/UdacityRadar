clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8

frequencyOfOperation = 77e9;
maxRange = 200;
rangeResolution = 1;
maxVelocity = 100;

c = 3*10^8;

%% User Defined Range and Velocity of target
 
% velocity between -70 m/s and +70 m/s
% distance between 0 and maxRange
targetInitialDistance = 100; 
targetInitialVelocity = 27; % m/s


%% FMCW Waveform Generation

bSweepHz = c / 2 * rangeResolution;
bSweepMHz = bSweepHz / 10^6;

disp("bSweepHz: " + bSweepHz + " Hz");

Tchirp = 5.5 * 2 * (maxRange / c);
TchirpMicro = Tchirp * 10^6;

disp("TChirp: " + TchirpMicro + " microseconds");

slope = bSweepHz / Tchirp;

disp("slope: " + slope/10e9 + " GHz/s");

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)             
    rangeOfTarget = targetInitialDistance + t(i)*targetInitialVelocity;

    % Tx(i) = cos(2 * PI * (fc * t + alpha * t^2 / 2))
    Tx(i) = cos(2 * pi * (fc * t(i) + slope * t(i)^2 / 2));

    % Rx (i) = cos(2 * pi * (fc * (t-tau) + alpha * (t-tau)^2 / 2 ))
    tau = rangeOfTarget * 2 /c;
    Rx(i) = cos(2 * pi * (fc * (t(i) - tau) + slope * (t(i) - tau)^2 / 2 ));
    
    Mix(i) = Tx(i) .* Rx(i);
end

%% RANGE MEASUREMENT

signal = reshape(Mix, [Nr Nd]);

signalFFT = fft(signal);

signalFFT = abs(signalFFT / Nr);

signalFFT = signalFFT(1:Nr/2+1);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

plot(signalFFT)

axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

Tr = 10;
Td = 8;

Gr = 4;
Gd = 4;

offset = 6;

thresholdBlock = zeros(size(RDM));

max_T = 1;
for i = (Tr + Gr) + 1 : Nr/2 - (Gr + Tr)
    for j = (Td + Gd) + 1 : Nd - (Gd + Td)
        noise_level = 0;

        for p = i - (Tr + Gr) : i + (Tr + Gr)
            for q = j - (Td + Gd) : j + (Td + Gd)
                if (abs(i-p) > Gr || abs(j-q) > Gd)
                    noise_level = noise_level + db2pow(RDM(p, q));
                end
            end
        end

        threshold = pow2db(noise_level / (2 * (Td + Gd + 1) * 2 * (Tr + Gr + 1) - (Gr * Gd) - 1));
        threshold = threshold + offset;
        
        CUT = RDM(i, j);
        if(CUT < threshold)
            thresholdBlock(i, j) = 0;
        else
            thresholdBlock(i, j) = max_T;
        end
    end
end

figure,surf(doppler_axis, range_axis, thresholdBlock);
colorbar;
