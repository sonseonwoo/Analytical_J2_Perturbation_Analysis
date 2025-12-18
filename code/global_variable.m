global Cnm Snm eopdata swdata SOLdata DTCdata APdata PC

load DE440Coeff.mat
PC = DE440Coeff;

% read Earth gravity field coefficients
Cnm = zeros(361,361);
Snm = zeros(361,361);
fid = fopen('GGM03C.txt','r');
for n=0:360
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end
fclose(fid);

% read Earth orientation parameters
fid = fopen('EOP-All.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    k = strfind(tline,'NUM_OBSERVED_POINTS');
    if (k ~= 0)
        numrecsobs = str2double(tline(21:end));
        tline = fgetl(fid);
        for i=1:numrecsobs
            eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        for i=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2double(tline(22:end));
        tline = fgetl(fid);
        for i=numrecsobs+1:numrecsobs+numrecspred
            eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        break
    end
end
fclose(fid);

% read space weather data
fid = fopen('SW-All.txt','r');
%  ---------------------------------------------------------------------------------------------------------------------------------
% |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
% | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
%  ---------------------------------------------------------------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    k = strfind(tline,'NUM_OBSERVED_POINTS');
    if (k ~= 0) 
        numrecsobs = str2double(tline(21:end));
        tline = fgetl(fid);
        for i=1:numrecsobs
            swdata(:,i) = fscanf(fid,'%i %d %d %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %f %i %i %f %i %f %f %f %f %f',[33 1]);
        end
        % remove the row of the Q parameter
        swdata = [swdata(1:27,:);swdata(29:33,:)];
        for i=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2double(tline(28:end));
        tline = fgetl(fid);
        %  -------------------------------------------------------------------------------------------------------------------------------
        % |                                                                                             Adj   Adj   Adj   Obs   Obs   Obs 
        % | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Ctr81 Lst81 F10.7 Ctr81 Lst81
        %  -------------------------------------------------------------------------------------------------------------------------------
        for i=numrecsobs+1:numrecsobs+numrecspred
            swdata(:,i) = fscanf(fid,'%i %d %d %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %f %i %i %f %f %f %f %f %f',[32 1]);
        end
        break
    end
end
fclose(fid);

% read solar storm indices
fid = fopen('SOLFSMY.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% read Ap data
fid = fopen('SOLRESAP.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD  F10 F10B Ap1 to Ap8
%  ------------------------------------------------------------------------
APdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d',[12 inf]);
fclose(fid);

% read geomagnetic storm indices
fid = fopen('DTCFILE.txt','r');
%  ------------------------------------------------------------------------
% | DTC YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
fclose(fid);
