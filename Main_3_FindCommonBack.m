clear all, close all
st = tic;
wavelength = 624; %nm
froot = 'Data/Cells/ProcessedData\MDAMB231\';
fstart = 'QPM10x_MDAMB231_';


numPos_Phase = [1];
numFrames_Phase = [7:16];


Btotal            = zeros(1200,1920);
c= zeros(1200,1920);
% c =0;
for pp = numPos_Phase
    tic
for ff = numFrames_Phase
    fname = [froot,fstart,'pos',num2str(pp,'%03d'),'_frame',num2str(ff,'%03d'),'.mat'];
    disp(fname)
        if isfile(fname)

            [D1,L1] = LoadSegment_short(fname, wavelength);
            BWfinal = L1 > 0;
            Btotal(~BWfinal) = (Btotal(~BWfinal).*c(~BWfinal) + D1(~BWfinal))./(c(~BWfinal)+1); 
            c(~BWfinal)=c(~BWfinal)+1;
        end

end
display(['Finished pos ',num2str(pp)])
toc
end

dim = 3;
save([froot, 'Btotal'], 'Btotal');
toc

%%

clear all, close all
st = tic;
wavelength = 624; %nm
froot = 'ProcessedData\MTG021_MTG084\';
fstart = 'QPM10x_MTG_';


numPos_Phase = [3];
numFrames_Phase = [5:14];


Btotal            = zeros(1200,1920);
c= zeros(1200,1920);
% c =0;
for pp = numPos_Phase
    tic
for ff = numFrames_Phase
    fname = [froot,fstart,'pos',num2str(pp,'%03d'),'_frame',num2str(ff,'%03d'),'.mat'];
    disp(fname)
        if isfile(fname)

            [D1,L1] = LoadSegment_short(fname, wavelength);
            BWfinal = L1 > 0;
            Btotal(~BWfinal) = (Btotal(~BWfinal).*c(~BWfinal) + D1(~BWfinal))./(c(~BWfinal)+1); 
            c(~BWfinal)=c(~BWfinal)+1;
        end

end
display(['Finished pos ',num2str(pp)])
toc
end

dim = 3;
save([froot, 'Btotal'], 'Btotal');
toc
