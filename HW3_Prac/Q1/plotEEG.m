function  plotEEG(X,text)
%X  % 8-channel data
% Plot Data
% Use function disp_eeg(X,offset,feq,ElecName)
offset = max(abs(X(:))) ;
feq = 100 ;
ElecName = [];
disp_eeg(X,offset,feq,ElecName,text);
end