function [phases] = hilbert_phases(Dat)

% Computation of instantaneous phases from Hilbert Transformation
% 10% of data length eliminated in both side in the Chris' code.

% L = size(Dat,1);
hilb = hilbert(Dat); % giving the Hilbert transform of the data, which is complex numbers.
phases = unwrap(angle(hilb)); %phases = unwrap(angle(hilb(round(L/10):end-round(L/10),:)));
end