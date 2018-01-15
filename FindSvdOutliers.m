function [tolerance,thresh_entry] = FindSvdOutliers(LHS,string)
% FINDSVDOUTLIERS search for the singular values (SV) of the LHS
% that are detrimental for the solution of the system, even if they do not
% cause the reciprocal of the condition number to be very low. 
% 
% FindSvdOutliers is called by MAIN**
% Input:
% the LHS of the solving system and "string", that defines the algortihm 
% to be applied:
%  if string = 'fit' -  a linear fit of the (logarithm of the) first 
%    'rng'*100% of the singular values of the LHS matrix is computed.
%    The standard deviation from this linear fit is computed next. Then,
%    the routine looks for singular values that are far (more than 'ndev' 
%    times the standard deviation) from the linear fit and gets the 
%    logarithm of the first singular value where this happens. It returns 
%    this value rounded towards the ceil to be used for setting the pseudo-
%    inverse tolerance.
%  if string = 'jump' -  the differences between successive values 
%    of the (logarithm of the) first 'rng'*100% of the singular values of 
%    the LHS matrix are computed and then the absolute maximum of these 
%    differences are taken. The routine looks in the remaining range of 
%    log(SVD) data and checks if there are differences larger than 
%    'njumps' * (the maximum value of the first 'rng'*100% singular 
%    values). It then gets the logarithm of the first singular value where 
%    this happens. It returns this value rounded towards the ceil to be 
%    used for setting the pseudo-inverse tolerance.
% 
% Output/Returns to MAIN***:
%  thresh_entry the first "outlier" sigular value and tolerance to be used 
%  as a singular value truncation tolerance
%
% BIBLIOGRAPHY
% 1. FreeHyTE Page - https://sites.google.com/site/ionutdmoldovan/freehyte
% 2. Moldovan ID, Cismasiu I - FreeHyTE: theoretical bases and developer’s 
% manual, https://drive.google.com/file/d/0BxuR3pKS2hNHTzB2N2Q4cXZKcGc/view
% 3. FreeHyTE Heat HTTE User's Manual - 
%    https://drive.google.com/drive/folders/0BxuR3pKS2hNHaFhiSjZHOE9TMzg
% 4. Geraldes, MRCB,  Elementos finitos híbridos-Trefftz adaptativos para
%   problemas de condução de calor. MSc Thesis, Universidade Nova de Lisboa,
%   2016 (in Portuguese).
%
% The outlier singular value search procedure implemented here is
% documented in reference 2 (Section 6.2.2) and reference [4] (Section
% 4.4.4).

%% The 'fit' procedure
if strcmp(string,'fit')
    
    % Basically, the 'fit' option computes a linear fit of the (logarithm
    % of the) first 'rng'*100% of the singular values of the LHS matrix
    % and then it computes the standard deviation from this linear fit on the
    % same data. Then, the routine looks for singular values that are
    % far (more than 'ndev' times the standard deviation) from the linear
    % fit and gets the logarithm of the first singular value where this
    % happens. It returns this value rounded towards the ceil to be used
    % for setting the pseudo-inverse tolerance.
    
    %% Initialization
    
    % The number of standard deviations to call an outlier (ndev) and the
    % 'healthy' data range for fitting (rng) are established empirically.
    % Since unstable solutions may compromise the convergence of the
    % iterative process, the parameters are highly biased towards accepting
    % false positives (identification of possibly faulty solving systems
    % when this is not, in fact, the case) than false negatives (failures
    % to identify system instability).
    ndev = 8;    
    rng = 0.85;  
    
    %% Statistical SV analysis
    % Computing the singular values of the LHS matrix...
    SV=svd(LHS); 
    % ...and their logarithms
    L10SV=log10(SV); 
    
    % index vector of the SV
    ind=1:length(L10SV); 
    ind=ind';
    
    % data range considered healthy
    okdata = round(rng*length(ind)); 
    
    % Computes linear fit on the healthy region okdata
    p=polyfit(ind(1:okdata),L10SV(1:okdata),1); 
    fitval = polyval(p,ind);
    
    % Computes difference between exact and linearly fitted data on the 
    % healthy range...
    dif = fitval(1:okdata)-L10SV(1:okdata); 
    
    % ...and standard deviation
    deviation = std(dif); 
    
    % Computes the deviation for a SV to be considered an outlier
    outval = ndev*deviation; 
    
    %% SV outlier identification
    % Computes difference between exact and linearly fitted data on the 
    % whole range of SV
    dif = fitval-L10SV;
    
    % identify the first and the last entries in dif that lay above outval. 
    % Note that dif is not used with its absolute value in order to only 
    % identify outliers that are below the linear fit. This means that we
    % are actually looking only for DROPS in the SV alignment.
    thresh_entry =  find(dif > outval,1,'first');
    thresh_last = find(dif > outval,1,'last');
    
    % Identifies the first and last singular values in the outlier region
    thresh_val = L10SV(thresh_entry);
    last_val = L10SV(thresh_last);
    
    % If the smallest outlier SV is larger than 1.E-6, no red flag is risen
    % by the routine. The reason for this is that for very small solving
    % systems, the range of 'healthy' data is insufficient for its
    % statistical analysis to be relevant. SV outliers are therefore not
    % identified if the 'worst' of them is not smaller than 1.E-6.
    if last_val < -6
        % the tolerance for trimming the SV is rounded towards the ceiling
        tolerance = ceil(thresh_val*10)/10; 
    else
        tolerance = [];
    end
    
    %% Debug mode - plots the SV and identifies the outlier
%     if ~isempty(tolerance) 
%         Fig=figure; plot(L10SV,'k*');hold on;%plot(fitval,'b');plot(dif,'r');
%         plot(thresh_entry, L10SV(thresh_entry),'o','markers',12); grid on;
%     end
%     set(Fig,'name',sprintf('Iteration %i',iteration));
%     set(Fig,'numbertitle','off');
    



%% The 'jump' procedure
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
    
    %% Initialization
    
    % The maximum jump multiplier to call an outlier (njumps) and the
    % 'healthy' data range for fitting (rng) are established empirically.
    % Since unstable solutions may compromise the convergence of the
    % iterative process, the parameters are highly biased towards accepting
    % false positives (identification of possibly faulty solving systems
    % when this is not, in fact, the case) than false negatives (failures
    % to identify system instability).
    njumps = 1.2;    
    rng = 0.9;  
    
    %% SV analysis
    % Computing the singular values of the LHS matrix...
    SV=svd(LHS); 
    % ...and their logarithms
    L10SV=log10(SV); 
    
    % total number of SVD
    nosvd=length(L10SV); 
    
    % data range considered healthy
    okdata = round(rng*nosvd); 
    
    % Computing the differences between consecutive SV in the health range. 
    % Since the SV are listed from the largest to the smallest, the gaps 
    % are negative. MIN is used, therefore, to find their maximum absolute
    % value.
    outval = min(L10SV(2:okdata+1)-L10SV(1:okdata)); 
    
    % Computing the differences between consecutive SV over the whole SV
    % range.
    dif = L10SV(okdata+1:end)-L10SV(okdata:end-1);
    
    % Identify the first and the last entries in dif that lay above outval. 
    % Note that dif is not used with its absolute value in order to only 
    % identify outliers that are below the linear fit. This means that we
    % are actually looking only for DROPS in the SV alignment.
    thresh_entry =  find(dif < njumps*outval,1,'first')+okdata;
    thresh_last = find(dif > outval,1,'last')+okdata;
    
    % Identifies the first and last singular values in the outlier region
    thresh_val = L10SV(thresh_entry);
    last_val = L10SV(thresh_last);
    
    % If the smallest outlier SV is larger than 1.E-6, no red flag is risen
    % by the routine. The reason for this is that for very small solving
    % systems, the range of 'healthy' data is insufficient for its
    % statistical analysis to be relevant. SV outliers are therefore not
    % identified if the 'worst' of them is not smaller than 1.E-6.    
    if last_val < -6
        % the tolerance for trimming the SV is rounded towards the ceiling        
        tolerance = ceil(thresh_val*10)/10; 
    else
        tolerance = [];
    end    
    
    %% Debug mode - plots the SV and identifies the outlier
%     if ~isempty(tolerance)
%         Fig=figure; plot(L10SV,'k*');hold on;plot(L10SV(2:end)-L10SV(1:end-1),'r');
%         plot(thresh_entry, L10SV(thresh_entry),'o','markers',12); grid on;
%     end
%     set(Fig,'name',sprintf('Iteration %i',iteration));
%     set(Fig,'numbertitle','off');
    
end
