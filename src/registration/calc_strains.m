function varargout = calc_strains(transform_file,filename_EV1,results_path,varargin)
%CALC_STRAINS calculates the full spatial jacobian from an elastix
% transform_file and then calculates maps of along fibre stress, along
% fibre shear and cross fibre_shear using the primary eigenvector map. The
% results are saved as NIfTI files.
%
% Bart Bolsterlee, Neuroscience Research Australia
% December 2017
%
% This function assumes that the command 'transformix' is recognised as an
% external command.
%
% USAGE: calc_strains(transform_file,EV1,results_path)
% or     calc_strains(transform_file,EV1,results_path,'mask',<mask_file>)
%
% INPUT:
% - transform_file : filename of the transformation file generated by
%                    Elastix.
% - filename_EV1   : filename of the NIfTI file with fibre orientation data
%                    (i.e. the primary eigenvector map as created with FSL
%                     or with DSI Studio using MakeEV1_map.m). All
%                     strain/deformation maps will be calculated at the
%                     resolution of the EV1 map.
% - results_path   : folder to which the resulting maps of
%                    strains/deformation/shear etc. are saved to.
%
% Optional inputs, provided as 'argument',<value> pairs:
% - mask     : filename of a binary mask. All voxels with value zero in the
%              mask will be set to zero in the resulting images.
%
% Note:  The data is presented in the NIfTI coordinate system, with x and y
%        axes flipped relative to the ITK coordinates in which Elastix
%        works.
%
% The following files will be saved to the results folder:
%
% - LAMBDA.nii.gz : Along fibre stress (lambda in eq. 3 in Blemker 2005)
% - B1.nii.gz : Along fibre shear (B1 in eq. 3 in Blemker 2005)
% - B2.nii.gz : Cross fibre shear (B2 in eq. 3 in Blemker 2005)
% - D.nii.gz  : length component of the shear (?)
% - deformation.nii.gz : deformation map
% - PSdir.nii.gz : map of principal strain direction
% - PSmag.nii.gz : map of principal strain magnitude
%
% The symbols and definitions used in this script are chosen to coincide
% with the descriptions on the following websites and book:
%
% http://www.continuummechanics.org/deformationgradient.html
% http://www.continuummechanics.org/polardecomposition.html
%
% Rheology: principles, measurements, and applications. Christopher W.
% Macosko, Wiley-VCH, 1994.
%
% The definition of the strain invariants and calculation of along-fibre
% stretch, along-fibre shear and cross-fibre shear can be found here:
%
% Blemker S, Pinsky P, Delp S (2005) A 3D model of muscle reveals the causes
% of nonuniform strains in the biceps brachii. J Biomech, 38.
%
% Criscione JC, Douglas AS, Hunter WC (2001) Physically based strain
% invariant set for materials exhibiting transversely isotropic behavior. J
% Mech Phys Solids, 49, 871-897.
%
% Pamuk U, Karakuzu A, Ozturk C, Acar B, Yucesoy CA (2016) Combined
% magnetic resonance and diffusion tensor imaging analyses provide a
% powerful tool for in vivo assessment of deformation along human muscle
% fibers. J Mech Behav Biomed Mater, 63, 207-19.

p = inputParser;
addRequired(p,'transform_file',@(x) exist(x,'file') == 2)
addRequired(p,'filename_EV1',@(x) contains(x,'.nii'))
addRequired(p,'results_path')
addParameter(p,'mask',[],@(x) contains(x,'.nii'))
addParameter(p,'F',[])
parse(p,transform_file,filename_EV1,results_path,varargin{:});

filename.mask    = p.Results.mask;
F               = p.Results.F;

% Create temporary working directory.
char_list = char(['a':'z' '0':'9']) ;
tmpdir = [];
while exist(tmpdir,'dir') == 7 || isempty(tmpdir)
    tmpdir = fullfile(tempdir,char_list(ceil(length(char_list)*rand(1,8))));
end
mkdir(tmpdir)

% Load the mask, if provided.
if ~isempty(filename.mask)
    mask = load_untouch_nii(filename.mask);
else
    mask = [];
end

try
    % Load the EV1 map, then change sign of x and y direction. I've checked
    % this in ITK-snap but may need some more checking (especially for
    % different srow combinations).
    
    EV1 = load_untouch_nii(filename_EV1);
%     EV1.img(:,:,:,1:2) = -EV1.img(:,:,:,1:2);
    

    if isempty(F)
        % Copy the transform file to the local working directory and modify the
        % spatial information to match the reference image (if a reference
        % images is provided).
        transform_tmp = fullfile(tmpdir,'transform.txt');
        if ~isempty(filename_EV1)
            ModifyTransformFile(transform_file,EV1,transform_tmp)
        else
            copyfile(transform_file,transform_tmp)
        end

        % Calculate deformation field
        transformix_cmd = sprintf('transformix -def all -out %s -tp %s',...
            tmpdir,transform_tmp);
        system(transformix_cmd);

        def = load_untouch_nii(fullfile(tmpdir,'deformationField.nii.gz'));

        % Flip coordinates to NIfTI coordinates
    %     def.img(:,:,:,:,1:2) = -def.img(:,:,:,:,1:2);

        % Calculate full spatial jacobian
        transformix_cmd = sprintf('transformix -jacmat all -out %s -tp %s',...
            tmpdir,transform_tmp);
        system(transformix_cmd);

        % Load spatial jacobian.
        F     = load_untouch_nii(fullfile(tmpdir,'fullSpatialJacobian.nii.gz'));
    end
    % Again, flip to NIfTI coordinates.
    % Elastix works in ITK coordinates, which have x and y-axis in opposite
    % directions to the NIfTI coordinates in which we work here. Correct
    % for that here by changing the signs of some of the componenets of the
    % Jacobian matrix (only the components that have dz in it (except dz/dz)
    % because for dx/dx and dy/dz for example the negative signs cancel out).
%     F.img(:,:,:,:,[3 6 7 8]) = -F.img(:,:,:,:,[3 6 7 8]);
    
    % Get image dimensions
    imdim = EV1.hdr.dime.dim(2:4);
    
    %% Calculate strain tensors
    
    % This needs to checked. It could also be that F21 = F(:,:,:,4) instead of
    % F12 = F(:,:,:,4), etc. This will depend on how the values of the jacobian
    % matrix are stored in the NIfTI file.
    %
    %     [F11 F12 F13]   [1 4 7]   [dx/dx dx/dy dx/dz]
    % F = [F21 F22 F23] = [2 5 8] = [dy/dx dy/dy dy/dz]
    %     [F31 F32 F33]   [3 6 9]   [dz/dx dz/dy dz/dz]
    %
%     F11 = F.img(:, :, :, 1,1); F12 = F.img(:, :, :, 1,4); F13 = F.img(:, :, :, 1,7);
%     F21 = F.img(:, :, :, 1,2); F22 = F.img(:, :, :, 1,5); F23 = F.img(:, :, :, 1,8);
%     F31 = F.img(:, :, :, 1,3); F32 = F.img(:, :, :, 1,6); F33 = F.img(:, :, :, 1,9);
%     
%     % Remove dilatoric component by dividing by the determinant of F.
%     
%     
%     % Calculate the Right Cauchy-Green Deformation Tensor as C = F'*F
%     %     [1 4 7]
%     % C = [2 5 8]
%     %     [3 6 9]
%     C = zeros([imdim 3 3]);
%     C(:, :, :, 1,1) = F11 .* F11 + F21 .* F21 + F31 .* F31;
%     C(:, :, :, 2,1) = F12 .* F11 + F22 .* F21 + F32 .* F31;
%     C(:, :, :, 3,1) = F13 .* F11 + F23 .* F21 + F33 .* F31;
%     C(:, :, :, 1,2) = F11 .* F12 + F21 .* F22 + F31 .* F32;
%     C(:, :, :, 2,2) = F12 .* F12 + F22 .* F22 + F32 .* F32;
%     C(:, :, :, 3,2) = F13 .* F12 + F23 .* F22 + F33 .* F32;
%     C(:, :, :, 1,3) = F11 .* F13 + F21 .* F23 + F31 .* F33;
%     C(:, :, :, 2,3) = F12 .* F13 + F22 .* F23 + F32 .* F33;
%     C(:, :, :, 3,3) = F13 .* F13 + F23 .* F23 + F33 .* F33;
%     clear F11 F12 F13 F21 F22 F23 F31 F32 F33
    
    %     % Or the Left Cauchy-Green Deformation Tensor : B = F * F'
    %     %     [1 4 7]
    %     % B = [2 5 8]
    %     %     [3 6 9]
    %     B = zeros([imdim 3 3]);
    %     B(:, :, :, 1,1) = F11 .* F11  +  F12 .* F12  +  F13 .* F13;
    %     B(:, :, :, 2,1) = F21 .* F11  +  F22 .* F12  +  F23 .* F13;
    %     B(:, :, :, 3,1) = F31 .* F11  +  F32 .* F12  +  F33 .* F13;
    %     B(:, :, :, 1,2) = F11 .* F21  +  F12 .* F22  +  F13 .* F23;
    %     B(:, :, :, 2,2) = F21 .* F21  +  F22 .* F22  +  F23 .* F23;
    %     B(:, :, :, 3,2) = F31 .* F21  +  F32 .* F22  +  F33 .* F23;
    %     B(:, :, :, 1,3) = F11 .* F31  +  F12 .* F32  +  F13 .* F33;
    %     B(:, :, :, 2,3) = F21 .* F31  +  F22 .* F32  +  F23 .* F33;
    %     B(:, :, :, 3,3) = F31 .* F31  +  F32 .* F32  +  F33 .* F33;
    %     clear F11 F12 F13 F21 F22 F23 F31 F32 F33
    
    
    %% Calculate strain invariants
    %
    hwait = waitbar(0,'Calculating along-fibre stretch, along-fibre shear and cross-fibre shear...',...
        'Name','Progress bar calc_strains');
    
    % Make empty maps for along fibre strain (LAMBDA), along-fibre shear (B1)
    % and cross-fibre shear (B2).
    LAMBDA = zeros(imdim);
    B1     = zeros(imdim);
    B2     = zeros(imdim);
    PSmag  = zeros(imdim);
    PSdir  = zeros([imdim,3]);
    detF   = zeros(imdim);
    for i = 1 : imdim(1)
        waitbar(i / imdim(1),hwait);
        for j = 1 : imdim(2)
            for k = 1 : imdim(3)
                
                if ~isempty(mask)
                    if mask.img(i,j,k) == 0
                        continue
                    end
                end
                % Get strain tensor and fibre direction of the current voxel
%                 Cvoxel = squeeze(C(i,j,k,1:3,1:3)) - eye(3);
                fib_dir = squeeze(EV1.img(i,j,k,:));
                if all(fib_dir==0);continue;end
                
                % Calculate the stretch tensor
%                 U = sqrtm(Cvoxel+eye(3));
                
                if isstruct(F)
                    % Get the deformation gradient (spatial jacobian)
                    Fvoxel = reshape(squeeze(F.img(i,j,k,:,:)),3,3);
                else
                    Fvoxel = F;
                end
                
                detF(i,j,k) = det(Fvoxel);
                % Remove the dilatoric component by dividing J^(1/3), where
                % J is the determinant of the deformation gradient matrix.
%                 Fvoxel = Fvoxel * det(Fvoxel)^(-1/3);
                
                % Calculate the Right Cauchy-Green Deformation Tensor as C = F'*F
                Cvoxel = Fvoxel' * Fvoxel;
                
                % Calculate principal strain magnitudes and directions
                % First, do an eigen-value decomposition on C.
                [EVec, EVal] = eig(Cvoxel);
                
                % Get the eigenvector associated with the largest
                % eigenvalue
                [max_EVal,max_idx] = max(diag(EVal));
                PSmag(i,j,k)       = max_EVal;
                PSdir(i,j,k,1:3)   = EVec(:,max_idx);
                
                % Calculate strain invariants. Definitions according to
                % Criscione et al. 2001, eq. 1.3.
                I1 = trace(Cvoxel); % trace of C
                I2 = (trace(Cvoxel)^2 - trace(Cvoxel^2))/2;
                I3 = det(Cvoxel);
                I4 = fib_dir' * Cvoxel * fib_dir;
                I5 = fib_dir' * Cvoxel^2 * fib_dir;
                
%                 % Definition in Pamuk et al. 2016, J Mech Beh Biom Mat.
%                 % This is indeed equivalent to definition by Blemker.
%                 J = det(Cvoxel);
%                 Lm = J^(-1/3)*sqrt((fib_dir'*Cvoxel*fib_dir));
%                 psi = sqrt((fib_dir'*Cvoxel^2*fib_dir) / (fib_dir'*Cvoxel*fib_dir)^2-1);                
%                 I4y = J^(2/3)*Lm^2;
%                 I5y = J^(4/3)*Lm^4*(1+psi^2);
                
                % Calculate along fibre stretch. Definition according to Blemker
                % et al. 2005, eq. 3.
                LAMBDA(i,j,k) = sqrt(I4);
                B1(i,j,k) = sqrt(I5 /(I4.^2)-1);
                % definition according to Blemker often gives imaginary
                % resutls.
%                 B2(i,j,k) = acosh( (I1*I4-I5) / (2*sqrt(I4)));
                
                % This is the definition from Criscione (2001) (variable
                % beta3 in eq. 5.3c).
                B2(i,j,k) = log( (I1*I4-I5) / (2*sqrt(I3*I4)) + sqrt( (I1*I4-I5)^2 / (2*sqrt(I3*I4))^2 - 1));
            end
        end
    end
    
    % Calculate the component of shear in the fibre direction (i.e. the
    % contribution of shear to lengthening)
%     D = B1 .* LAMBDA;
    close(hwait)
    
%     % Make sure that all principal strain vectors point in the positive z-direction
    PSdir = PSdir .* repmat(sign(PSdir(:,:,:,3)),1,1,1,3);
    %% Save the results
    
    % Save the along-fibre stress map as LAMBDA.nii.gz
    S = EV1;
    S.hdr.dime.dim = [3 imdim 1 1 1 1];
    S.img = cast(LAMBDA,'like',EV1.img);
    S.hdr.dime.glmax = max(S.img(:));
    S.hdr.dime.glmin = min(S.img(:));
    if exist(results_path,'dir') ~= 7
        mkdir(results_path)
    end
    filename.lambda = fullfile(results_path,'LAMBDA.nii.gz');
    save_untouch_nii(S,filename.lambda)
    fprintf('Along-fibre stretch field saved as %s\n',filename.lambda)
    
    % Save the along-fibre shear map as B1.nii.gz
    S.img = cast(B1,'like',EV1.img);
    S.hdr.dime.glmax = max(S.img(:));
    S.hdr.dime.glmin = min(S.img(:));
    filename.B1 = fullfile(results_path,'B1.nii.gz');
    save_untouch_nii(S,filename.B1)
    fprintf('Along-fibre shear field saved as %s\n',filename.B1)
    
    % Save the cross-fibre shear map as B2.nii.gz
    S.img = cast(B2,'like',EV1.img);
    S.hdr.dime.glmax = max(S.img(:));
    S.hdr.dime.glmin = min(S.img(:));
    filename.B2 = fullfile(results_path,'B2.nii.gz');
    save_untouch_nii(S,filename.B2)
    fprintf('Cross-fibre shear field saved as %s\n',filename.B2)
    
    % Save the determinant of F as detF.nii.gz
    S.img = cast(detF,'like',EV1.img);
    S.hdr.dime.glmax = max(S.img(:));
    S.hdr.dime.glmin = min(S.img(:));
    filename.detF = fullfile(results_path,'detF.nii.gz');
    save_untouch_nii(S,filename.detF)
    fprintf('Determinant of F saved as %s\n',filename.detF)
    
%     % Save the shear-length-component map as D.nii.gz
%     S.img = cast(D,'like',EV1.img);
%     S.hdr.dime.glmax = max(S.img(:));
%     S.hdr.dime.glmin = min(S.img(:));
%     filename.D = fullfile(results_path,'D.nii.gz');
%     save_untouch_nii(S,filename.D)
%     fprintf('Shear-length component field saved as %s\n',filename.D)
    
    % Save the principal strain magnitude as PSmag.nii.gz
    S.img = cast(PSmag,'like',EV1.img);
    S.hdr.dime.glmax = max(S.img(:));
    S.hdr.dime.glmin = min(S.img(:));
    filename.PSmag = fullfile(results_path,'PSmag.nii.gz');    
    save_untouch_nii(S,filename.PSmag)
    fprintf('Principal strain magnitude map saved as %s\n',filename.PSmag)
    
    % Save the principal strain direction as PSdir.nii.gz
    S.img = cast(PSdir,'like',EV1.img);
    S.hdr.dime.dim(1) = 4;
    S.hdr.dime.dim(5) = 3;
    S.hdr.dime.glmax = max(S.img(:));
    S.hdr.dime.glmin = min(S.img(:));
    filename.PSdir = fullfile(results_path,'PSdir.nii.gz');    
    save_untouch_nii(S,filename.PSdir)
    fprintf('Principal strain direction map saved as %s\n',filename.PSdir)
    
    if exist('def','var') == 1
        % Save the deformation field
        if ~isempty(mask)
            % Set values outside mask to zero.
            def.img = def.img .* repmat(cast(mask.img,'like',def.img)~=0,1,1,1,1,3);
        end
        filename.def = fullfile(results_path,'deformation.nii.gz');
        save_untouch_nii(def,filename.def)
        fprintf('Deformation field saved as as: %s\n',filename.def)
    end
    
    % Remove temporary working directory.
    rmdir(tmpdir,'s')
    
    % Return filenames as ouput, if an output argument is provided.
    if nargout > 0
        varargout{1} = filename;
    end
    
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    error(ME.message)
end

end % of function