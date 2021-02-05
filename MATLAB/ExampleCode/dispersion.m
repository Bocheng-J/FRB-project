% dispersion() - First-mode dispersion correction of a finite  
% arbitrary signal in a cylindrical bar
%
%
% OPERATION:
% - Finds FFT of the signal
% - Corrects phase velocity and amplitude of each frequency
%   using method described by Tyas & Pope (2005)
% - Reconstructs signal using IFFT
% - (Frequencies above fa/c0 = 0.2619 stripped (d/L = 0.6) 
%   due to limitations of M1 correction)
%
% INPUTS:
% x        Zero-padded strain signal in time domain 
%          (1xN numeric)
% fs       Sampling frequency, Hz
% a        Bar radius, m
% c0       One-dimensional wave velocity of the bar, m/s
% E        Young's modulus of the bar, GPa
% z        Distance to correct over, m 
%          (+ve in direction of propagation)
%
% OUTPUTS:
% xStrain  Dispersion-corrected strain signal
% xStress  Dispersion-corrected stress signal, MPa



function [xStrain xStress] = dispersion(x,fs,a,c0,E,z)

    % Input signal
    N = length(x); % Number of elements in signal
    dt = 1/fs; % Time step, s 
    t = 0:dt:dt*(N-1); % Time base, s
    f = (0:N-1)*(fs/N); % FFT frequencies, Hz
    fMax = 0.2619*c0/a; % Max correctable frequency due to factor M1 limitations, Hz

    % FFT the signal
    X = fft(x);
    XStrain = X; % Create a copy for strain correction
    XStress = X; % Create a copy for stress correction

    % Phase shift, adjust magnitude of frequency components
    numberOfBins = length(X);
    DCbin = 1; % DC frequency bin

    if(mod(numberOfBins,2)==0)
        % N is even
        positiveBins = 2:numberOfBins/2; % Positive frequency bins
        nyquistBin = numberOfBins/2+1; % Nyquist frequency bin
        binsToEdit = [positiveBins nyquistBin]; % Total bins to edit individually
        negativeBins = numberOfBins/2+2:numberOfBins; % Negative frequency bins (populated through conjugation of positive bins)
    else
        % N is odd
        positiveBins = 2:(numberOfBins+1)/2; % Positive frequency bins
        binsToEdit = positiveBins; % Total bins to edit individually
        negativeBins = (numberOfBins+1)/2+1:numberOfBins; % Negative frequency bins (populated through conjugation of positive bins)
    end

    for b = binsToEdit
        if f(b) <= fMax
            % Edit phase and amplitude of positive bins
            [angleMod M1 M2] = dispersionFactors(f(b),a,c0,z); % Find phase shift and factors M1 and M2 for current freq
            XStrain(b) = M1*abs(X(b)) * exp(1i * (angle(X(b))-angleMod)); % Apply phase shift and factors M1 to obtain corrected strain
            XStress(b) = M1*M2*abs(X(b)) * exp(1i * (angle(X(b))-angleMod)); % Apply phase shift and factors M1 and M2 to obtain corrected stress(/E)
        else
            % Above fMax zero X data (apply perfect low-pass filter)
            XStrain(b) = 0;
            XStress(b) = 0;
        end
    end

    XStrain(negativeBins) = conj(XStrain(positiveBins(end:-1:1))); % Correct negative bins by taking complex conjugate of positive bins
    XStress(negativeBins) = conj(XStress(positiveBins(end:-1:1))); % Correct negative bins by taking complex conjugate of positive bins

    % Convert the corrected frequency components back into the time domain
    xStrain = ifft(XStrain); % Corrected strain
    xStress = ifft(XStress)*E*1000; % Corrected stress, MPa

end

% dispersionFactors() - Calculate corrections to amplitude 
% and phase angle to account for dispersion at a particular 
% frequency
%
% REQUIRES:
% <dataTable>  A .mat file containing pre-calculated vectors 
%              for a particular Poisson's ratio
%              - normFreqs: normalised frequencies (f*a/c0)
%              - vRatios: normalised velocities (c/c0)
%              - M1: amplitude factor M1
%              - M2: normalised amplitude factor M2 (M2/E)
%
% INPUTS:
% f            Frequency, Hz
% a            Bar radius, m
% c0           One-dimensional wave velocity of the bar, m/s
% z            Distance to apply correction over, m
%
% OUTPUTS:
% angleMod  Phase angle correction, rad
% M1        Correction for variation in response across bar 
%           cross-section
% M2        Correction for variation in ratio of axial stress
%           and axial strain (dynamic Young's modulus)

function [angleMod M1 M2] = dispersionFactors(f,a,c0,z)
    load('dispersionFactors_PR29.mat'); % File containing phase velocity data
    normFreq = f*a/c0; % Normalised frequency

	% Find change in phase angle
	phaseVelocity = interp1(normFreqs,vRatios,normFreq)*c0; % Interpolated phase velocity value
    angleMod = 2 * pi * f * z / phaseVelocity; % Change in phase angle at normFreq
    
    % Find amplitude factors M1 and M2
	M1 = interp1(normFreqs,M1,normFreq); % Interpolated value of M1 at normFreq
	M2 = interp1(normFreqs,M2,normFreq); % Interpolated value of M2/E at normFreq
end