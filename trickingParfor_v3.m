% this function was originally created to dance around problems with
% putting the parfor loop in this file.

% version 2 converts the single precision phase data into a 16 bit integer
% to decrease the file size and shorten the save time, load time, and file
% transfer time.

function trickingParfor_v3(froot, savefolder, fstart, pos, frame)
%     F = @(x) fftshift(fft2(ifftshift(x)));
%     Ft = @(x) fftshift(ifft2(ifftshift(x)));
%         ii = 81;
%     pos = jj;
    % where you put the tiff files
    filedir = (froot);

    % setup output folder here
%     out_dir = (froot);
    %         Next line is commented because output directory should already
    %         exist
    %         mkdir(out_dir);

    % find all the figures
    imglist = dir([filedir,'Image_*.tiff']);
    %% images from directory    
    nfiles = length(imglist);    % Number of files found

    Ih = cell([nfiles,1]);

    for kk=1:nfiles
       currentfilename = [filedir,imglist(kk).name];
       %current = dir([filedir,currentfilename]);
       currentimage=double(imread(currentfilename));
       Ih{kk} = currentimage;
    end
    %         gets tiff for top image
    D = dir(currentfilename);
    timestamp = D.datenum;

    % load system parameters
    SystemSetupDPC();

    % what's the angles used in the DPC measurements?
            rot_an = [0 90];
%     rot_an = [90 180];

    nAngle = length(rot_an);

    %% calculating DPC Transfer function
    Hi = zeros(Np(1),Np(2),nAngle);

    % what's the illumination NA used?
    NA_illum = NA;
    source_led = 1;
    for dd = 1:2
        % source shape?
        S = Dsource_LR(rot_an(dd), NA_illum, lambda, u, v);
        S = S.*source_led;
        % calculate transfer function
        [~, Hi(:,:,dd)] = DPCTransFunc(pupil, S);
    end


       % calculate DPC image
    IDPC = zeros(Np(1),Np(2),nAngle);   
    for dd = 1:nAngle
        if dd==1
            I1 = Ih{((dd-1)*4)+3};
            I2 = Ih{((dd-1)*4)+2};
        end
        if dd>1
            I1 = Ih{((dd-2)*4)+1};
            I2 = Ih{((dd-2)*4)+4};
        end
        IDPC(:,:,dd) = (I1-I2)./(I1+I2);
    %             B=planefit(IDPC(:,:,dd));
    %             figure;
    %             imagesc(IDPC(:,:,dd));
    end
%     intensity = I1;


    % calculate DPC image
    % IDPC = zeros(Np(1),Np(2),nAngle);
    % for dd = 1:nAngle
    %     IDPC(:,:,dd) = (Ih(:,:,(dd-1)*2+1)-Ih(:,:,dd*2))./(Ih(:,:,(dd-1)*2+1)+Ih(:,:,dd*2));
    % end

    % NEED TO TUNE THE FOLLOWING REGULARIZATION PARAMETERS BASED ON
    % EXPERIMENTAL SNR
    reg2 = 4e-3;
    % calculate phase by deconvolution
    ph_dpc =  DPC_tik(IDPC, Hi, reg2);
    Phase = int16(ph_dpc(:,:)/(4*pi)*65536);
    save([savefolder, fstart,'pos',num2str(pos,'%03d'),'_frame',num2str(frame,'%03d'),], ...
    'Phase','timestamp')
%     imwrite(uint16(Phase),[savefolder, fstart,'pos',num2str(pos,'%03d'),'_frame',num2str(frame,'%03d'),'.tiff'],'tiff')
%     rmdir(filedir, 's');
