function [tolerance,thresh_entry] = FindSvdOutliers(LHS,string)
% The objective of this function is to find the singular values of the LHS
% that are detrimental for the solution of the system, even if they do not
% cause the reciprocal of the condition number to be very low.

if strcmp(string,'fit')
    
    % Basically, the 'fit' option computes a linear fit of the (logarithm
    % of the) first 'rng'*100% of the singular values of the LHS matrix
    % and then it computes the standard deviation from this linear fit on the
    % same data. Then, the routine looks for singular values that are
    % far (more than 'ndev' times the standard deviation) from the linear
    % fit and gets the logarithm of the first singular value where this
    % happens. It returns this value rounded towards the ceil to be used
    % for setting the pseudo-inverse tolerance.
    
    ndev = 8;    % the number of standard deviations to call an outlier
    rng = 0.85;  % the data range for fitting
    
    SV=svd(LHS); % singular values of the LHS matrix...
    L10SV=log10(SV); % ...and their logarithms
    ind=1:length(L10SV); % entry vector
    ind=ind';
    okdata = round(rng*length(ind)); % data range considered healthy
    p=polyfit(ind(1:okdata),L10SV(1:okdata),1); % linear fit on the healthy region
    fitval = polyval(p,ind);
    dif = fitval(1:okdata)-L10SV(1:okdata); % difference between exact and linearly fitted data on the healthy range
    deviation = std(dif); % standard deviation
    outval = ndev*deviation; % deviation to be considered an outlier
    dif = fitval-L10SV;
    % identify the first entry in dif which lays above outval. Note that dif is
    % not used with its absolute value to only identify thos outliers that are
    % below the linear fit.
    thresh_entry =  find(dif > outval,1,'first');
    thresh_last = find(dif > outval,1,'last');
    % the first singular value in the outlier region
    thresh_val = L10SV(thresh_entry);
    maxdif_val = L10SV(thresh_last);
    if maxdif_val < -6
        tolerance = ceil(thresh_val*10)/10; % round upwards to the first decimal
    else
        tolerance = [];
    end
    
%     if ~isempty(tolerance) 
%         Fig=figure; plot(L10SV,'k*');hold on;%plot(fitval,'b');plot(dif,'r');
%         plot(thresh_entry, L10SV(thresh_entry),'o','markers',12); grid on;
%     end
%     set(Fig,'name',sprintf('Iteration %i',iteration));
%     set(Fig,'numbertitle','off');
    
else
    
    % The 'jump' option computes the differences between successive values 
    % of the (logarithm of the) first 'rng'*100% of the singular values of 
    % the LHS matrix and then it takes the absolute maximum of these 
    % differences. Then, the routine looks in the remaining range of 
    % log(SVD) data and checks if there are differences larger than 
    % 'njumps' * (the maximum value so far). It then gets the logarithm of 
    % the first singular value where this happens. It returns this value 
    % rounded towards the ceil to be used for setting the pseudo-inverse 
    % tolerance.
    
    njumps = 1.2;    % the number of 'jumps' to call an outlier
    rng = 0.9;  % the data range for finding the maximum 'jump'
    
    SV=svd(LHS); % singular values of the LHS matrix...
    L10SV=log10(SV); % ...and their logarithms
    nosvd=length(L10SV); % total number of SVD
    okdata = round(rng*nosvd); % data range considered healthy
    outval = min(L10SV(2:okdata+1)-L10SV(1:okdata)); % min is used as dif values are negative
    dif = L10SV(okdata+1:end)-L10SV(okdata:end-1);
    
    % identify the first entry in dif which lays above outval. Note that dif is
    % not used with its absolute value to only identify thos outliers that are
    % below the linear fit.
    thresh_entry =  find(dif < njumps*outval,1,'first')+okdata;
    thresh_last = find(dif > outval,1,'last')+okdata;
    % the first singular value in the outlier region
    thresh_val = L10SV(thresh_entry);
    maxdif_val = L10SV(thresh_last);
    if maxdif_val < -6
        tolerance = ceil(thresh_val*10)/10; % round upwards to the first decimal
    else
        tolerance = [];
    end    
    
%     if ~isempty(tolerance)
%         Fig=figure; plot(L10SV,'k*');hold on;plot(L10SV(2:end)-L10SV(1:end-1),'r');
%         plot(thresh_entry, L10SV(thresh_entry),'o','markers',12); grid on;
%     end
%     set(Fig,'name',sprintf('Iteration %i',iteration));
%     set(Fig,'numbertitle','off');
    
end