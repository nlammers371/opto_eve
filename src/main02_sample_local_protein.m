% main02_sample_local_protein(project, rawPath, proteinChannel, varargin)
%
% DESCRIPTION
% Script to generate samples of local protein contentration at active gene
% loci and selected control locations
%
% ARGUMENTS
% project: master ID variable 
%
% rawPath: Full or relative path to PreProcessed folder (pr
%               equivalent)
%
% proteinChannel: Integer corresponding to protein channel
%
% OPTIONS
% dropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
% zeissFlag: Integer specifying a pixel correction (1 pixel needed for
%               data taken on Zeiss780)
%
% OUTPUT: nucleus_struct_protein: compiled data set with protein samples

function nucleus_struct_protein = main02_sample_local_protein(project,DropboxFolder,varargin)
addpath('./utilities')
ROIRadiusSpot = .2; % radus (um) of region used to query and compare TF concentrations
minSampleSepUm = 1.5; %um
mf_samp_rad = 0.8; % distance (um) from nucleus center to include in sample 
minEdgeSepUm = .5; %um
segmentNuclei = 0;
% PSF info for 3D sampling
use_psf_fit_dims = false; % if true, use fits from PSF fitting
xy_sigma_um = 0.25;% um 
z_sigma_um = 0.6; % um
ignoreQC = false;
zMax = 30;
rawPath = 'E:\LocalEnrichment\Data\PreProcessedData\';
proteinChannel = 1;
[~, DataPath, ~] =   header_function(DropboxFolder, project);

for i = 1:(numel(varargin)-1)  
    if i ~= numel(varargin)      
        eval([varargin{i} '=varargin{i+1};']);        
    end    
end

% Load trace data
load([DataPath '/nucleus_struct.mat'],'nucleus_struct')
if use_psf_fit_dims
    load([DataPath '/psf_dims.mat'],'psf_dims')
end
load([DataPath '/set_key.mat'],'set_key')
snipPath = [DataPath '/qc_images/'];
refPath = [DataPath '/refFrames/'];
mkdir(refPath)
mkdir(snipPath)

% get MCP channel
options = [1 2];
mcp_channel = options(options~=proteinChannel);
write_snip_flag = false;
%%%%%%%%%%%%%%%%%%%%%%%%%
% remove all frames that do not contain a segmented particle or that

threeD_flag = nucleus_struct(1).threeD_flag;
% threeD_flag = false;
% remove entries with no particle
nucleus_struct = nucleus_struct(~isnan([nucleus_struct.ParticleID]));

fnames = fieldnames(nucleus_struct);
for i = 1:numel(nucleus_struct)
    nucleus_struct(i).time_orig = nucleus_struct(i).time;
    fluo = nucleus_struct(i).fluo;
    nan_ft = ~isnan(fluo);
    for j = 1:numel(fnames)
        vec = nucleus_struct(i).(fnames{j});
        if numel(vec) == numel(fluo)        
            nucleus_struct(i).(fnames{j}) = vec(nan_ft);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% generate indexing vectors
frame_ref = [nucleus_struct.frames];
nc_x_ref = [nucleus_struct.xPos];
nc_y_ref = [nucleus_struct.yPos];

spot_x_ref = [nucleus_struct.xPosParticle];
spot_y_ref = [nucleus_struct.yPosParticle];
spot_z_ref = [nucleus_struct.zPosParticle];

if threeD_flag
    spot_x_ref3D = [nucleus_struct.xPosParticle3D];
    spot_y_ref3D = [nucleus_struct.yPosParticle3D];
    spot_z_ref3D = [nucleus_struct.zPosParticle3D];
else
    spot_x_ref3D = [nucleus_struct.xPosParticle];
    spot_y_ref3D = [nucleus_struct.yPosParticle];
    spot_z_ref3D = [nucleus_struct.zPosParticle];
end

set_ref = [];
master_ind_ref = [];
lin_ind_ref = [];
sub_ind_ref = [];
pt_ref = [];
pt_qc_ref = [];
for i = 1:numel(nucleus_struct)
    ParticleID = nucleus_struct(i).ParticleID;    
    pt_ref = [pt_ref repelem(ParticleID, numel(nucleus_struct(i).frames))];
    pt_qc_ref = [pt_qc_ref repelem(nucleus_struct(i).qc_flag, numel(nucleus_struct(i).frames))];
    set_ref = [set_ref repelem(nucleus_struct(i).setID,numel(nucleus_struct(i).frames))];
    master_ind_ref = [master_ind_ref repelem(round(nucleus_struct(i).ncID*1e5),numel(nucleus_struct(i).frames))];
    lin_ind_ref = [lin_ind_ref repelem(i,numel(nucleus_struct(i).frames))]; 
    sub_ind_ref = [sub_ind_ref 1:numel(nucleus_struct(i).frames)];
end

% option to override qc 
if ignoreQC 
    pt_qc_ref = true(size(pt_qc_ref));
end
%%%%%%%%%%%%%%%%%%%%
%%% Set size parameters
set_vec = [nucleus_struct.setID];
set_index = unique(set_vec);

PixelSize = nucleus_struct(1).PixelSize;
zStep = nucleus_struct(1).zStep;
% determine size of neighborhood to use
nb_size = round(10 ./ PixelSize);

% size of gaussian smoothing kernel 
sm_kernel = round(1 ./ PixelSize);

% set min and max acceptable area for nucleus segmentation
min_area = round(pi*(2 ./ PixelSize).^2);
max_area = round(pi*(4 ./ PixelSize).^2);

% set snippet to be 3um in size
pt_snippet_size = round(1.5 ./ PixelSize);

% set min separation between control and locus to 2um
minSampleSep = round(minSampleSepUm ./ PixelSize);
minEdgeSep = round(minEdgeSepUm ./ PixelSize);

% calculate ROI size in pixels for spot and control
roi_rad_spot_pix = round(ROIRadiusSpot ./ PixelSize);

% calculate average frame-over-frame particle drift from data
lin_diff_vec = diff(lin_ind_ref);
x_diff_vec = diff(spot_x_ref);
y_diff_vec = diff(spot_y_ref);
dr_vec = sqrt(x_diff_vec.^2+y_diff_vec.^2);
dr_vec = dr_vec(lin_diff_vec==0);
% sets sigma of movement for virtual spot
driftTol = nanmedian(dr_vec)*PixelSize;

% set dims for 3D protein sampling
if use_psf_fit_dims
    xy_sigma = psf_dims.xy_sigma;
    z_sigma = psf_dims.z_sigma;
else
    xy_sigma = round(xy_sigma_um/PixelSize,1);% pixels
    z_sigma = round(z_sigma_um/zStep,1); % pixels
end

%%%%%%%%%%%%%%%%%%%%%%%
% Designate fields ot be added to nucleus structure
new_vec_fields = {'spot_protein_vec_3d','spot_protein_vec', 'serial_null_protein_vec',...
    'serial_null_protein_vec_3d','edge_null_protein_vec','edge_null_protein_vec_3d','mf_null_protein_vec',...
    'spot_mcp_vec','edge_mcp_protein_vec','serial_qc_flag_vec','edge_qc_flag_vec', ...
    'edge_null_x_vec', 'serial_null_x_vec','serial_null_y_vec','edge_null_y_vec', 'edge_null_nc_vec',...
    'spot_edge_dist_vec','serial_null_edge_dist_vec'};


new_snip_fields = {'spot_protein_snips', 'edge_null_protein_snips',...
    'spot_mcp_snips','edge_null_mcp_snips'};

% Initialize fields
for i = 1:numel(nucleus_struct)
    ref = nucleus_struct(i).xPos;
    for j = 1:numel(new_vec_fields)
        nucleus_struct(i).(new_vec_fields{j}) = NaN(size(ref));
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%%% make source key
src_cell = {1,numel(set_index)};
% pt_id_vec = [nucleus_struct.ParticleID];
for s = 1:numel(set_index)
    ind = find([nucleus_struct.setID]==set_index(s),1);   
    src = nucleus_struct(ind).source_path;    
    src_cell{s} = src;     
end

%%% Generate reference array for set-frame combos
set_frame_array = unique([set_ref' frame_ref'],'row');
qc_structure = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Nucleus segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%
% first check to see if segmentation files exist
segment_indices = 1:size(set_frame_array,1);
spot_frame_vec = false(1,size(set_frame_array,1));
nc_frame_vec = false(1,size(set_frame_array,1));
for i = 1:size(set_frame_array,1)        
    setID_temp = set_frame_array(i,1);
    frame_temp = set_frame_array(i,2);  
    nc_ref_name = [refPath 'nc_ref_frame_set' sprintf('%02d',setID_temp) '_frame' sprintf('%03d',frame_temp) '.mat'];
    nc_frame_vec(i) = isfile(nc_ref_name);    
    spot_ref_name = [refPath 'spot_roi_frame_set' sprintf('%02d',setID_temp) '_frame' sprintf('%03d',frame_temp) '.mat'];    
    spot_frame_vec(i) = isfile(spot_ref_name);    
end    
if all(spot_frame_vec) && segmentNuclei   
    warning('previous segmentation results found')
    y = 1;
    n = 0;
    overwrite = input('overwrite segmentation results? (y/n)');
    segmentNuclei = overwrite;
elseif ~all(nc_frame_vec) && ~segmentNuclei   
    warning('some or all frames missing nucleus segmentation data. Segmenting missing frames only')
    segmentNuclei = 1;   
    segment_indices = find(~nc_frame_vec);
end


if segmentNuclei
    h = waitbar(0,'segmenting nuclei...');
%     disp('segmenting nuclei...')
    tic    
    %%% Segment nuclei     
    % initialize arrays to store segmentation info
    nucleus_frame_array = NaN(yDim,xDim,numel(segment_indices));
    spot_frame_array = NaN(yDim,xDim,numel(segment_indices));
    for w = 1:numel(segment_indices)
        waitbar(h,w/numel(segment_indices));
        i = segment_indices(w);
        setID_temp = set_frame_array(i,1);
        frame_temp = set_frame_array(i,2);  
        % get nucleus
        frame_set_filter = set_ref==setID_temp&frame_ref==frame_temp;
        nc_x_vec_temp = round(nc_x_ref(frame_set_filter));
        nc_y_vec_temp = round(nc_y_ref(frame_set_filter)); 
        % indexing vectors    
        nc_master_vec = master_ind_ref(frame_set_filter);  
        % get list of unique indices 
        [nc_master_vec_u,ia,~] = unique(nc_master_vec,'stable');
        % unique nucleus vectors
        nc_x_vec_u = round(nc_x_vec_temp(ia));
        nc_y_vec_u = round(nc_y_vec_temp(ia));    
        % particle positions        
        spot_x_vec = round(spot_x_ref(frame_set_filter));
        spot_y_vec = round(spot_y_ref(frame_set_filter));            

        src = set_key(set_key.setID==setID_temp,:).prefix{1};   
        protein_stack = load_stacks(rawPath, src, frame_temp, proteinChannel);
        % generate protein gradient frame for segmentation
        protein_smooth = imgaussfilt(mean(protein_stack,3),round(sm_kernel/2));                
        protein_grad = imgradient(protein_smooth);   
        % flatten background         
        protein_bkg = imgaussfilt(protein_smooth, round(nb_size/2));
        protein_grad_norm = protein_grad ./ protein_bkg;        
        
        % try local thresholding
        n_x_local = floor(xDim / nb_size);
        n_y_local = floor(yDim / nb_size);
        x_dim_local = round(xDim / n_x_local);
        y_dim_local = round(yDim / n_y_local);
        protein_bin_clean = false(size(protein_grad_norm));
        for x = 1:n_x_local
            for y = 1:n_y_local
                x_start = (x-1)*x_dim_local+1;
                x_stop = min([x*x_dim_local,xDim]);
                y_start = (y-1)*y_dim_local+1;
                y_stop = min([y*y_dim_local,yDim]);
                pt_section = protein_grad_norm(y_start:y_stop,x_start:x_stop);
                thresh = multithresh(pt_section);
                section_bin = pt_section > thresh; 
                section_bin_clean = bwareaopen(section_bin,sm_kernel^2);                
                % record
                protein_bin_clean(y_start:y_stop,x_start:x_stop) = bwmorph(section_bin_clean,'hbreak');
            end
        end        
        % label regions
        nc_frame_labels = logical(protein_bin_clean);         
        % frame info
        nc_lin_indices = sub2ind(size(protein_bin_clean),round(nc_y_vec_u),round(nc_x_vec_u));
        % take convex hull
        stats = regionprops(nc_frame_labels,'ConvexHull');
        nc_ref_frame = zeros(size(nc_frame_labels)); 
        for j = 2:numel(stats)
            hull_points = stats(j).ConvexHull;
            mask = poly2mask(hull_points(:,1),hull_points(:,2),yDim,xDim);   
            nc_bin_ids = mask(nc_lin_indices);
            if sum(nc_bin_ids) == 1 % enforce unique
                nc_ref_frame(mask) = nc_master_vec_u(nc_bin_ids);
            end
        end       
        % generate array indicating distance of each pixel from an active locus 
        nc_indices = sub2ind(size(nc_ref_frame),spot_y_vec,spot_x_vec);
        spot_dist_frame_temp = zeros(size(nc_ref_frame));
        spot_dist_frame_temp(nc_indices(~isnan(nc_indices))) = 1;
        spot_dist_frame_temp = bwdist(spot_dist_frame_temp);
        % label regions within designated integration radius of a spot    
        % store arrays
        nucleus_frame_array(:,:,w) = nc_ref_frame;
        spot_frame_array(:,:,w) = spot_dist_frame_temp;
        % save spot and nucleus reference frames    
    end
    disp('saving segmentation results...')
    % save arrays
    for w = 1:numel(segment_indices)
        i = segment_indices(w);
        setID_temp = set_frame_array(i,1);
        frame_temp = set_frame_array(i,2);  
        nc_ref_name = [refPath 'nc_ref_frame_set' sprintf('%02d',setID_temp) '_frame' sprintf('%03d',frame_temp) '.mat'];
        nc_ref_frame = nucleus_frame_array(:,:,w);
        save(nc_ref_name,'nc_ref_frame');
        spot_ref_name = [refPath 'spot_roi_frame_set' sprintf('%02d',setID_temp) '_frame' sprintf('%03d',frame_temp) '.mat'];
        spot_dist_frame = spot_frame_array(:,:,w);
        save(spot_ref_name,'spot_dist_frame');
    end
    disp('done.')
    toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Protein sampling
%%%%%%%%%%%%%%%%%%%%%%%%%
disp('taking protein samples...')
for i = 1:size(set_frame_array,1)    
    tic
    setID = set_frame_array(i,1);
    frame = set_frame_array(i,2);  
        
    % load spot and nucleus reference frames
    nc_ref_name = [refPath 'nc_ref_frame_set' sprintf('%02d',setID) '_frame' sprintf('%03d',frame) '.mat'];
    load(nc_ref_name,'nc_ref_frame');
    nc_dist_frame = bwdist(~nc_ref_frame);    
    spot_ref_name = [refPath 'spot_roi_frame_set' sprintf('%02d',setID) '_frame' sprintf('%03d',frame) '.mat'];
    load(spot_ref_name,'spot_dist_frame');
    spot_roi_frame = bwlabel(spot_dist_frame <= roi_rad_spot_pix);
    % get nucleus
    frame_set_filter = set_ref==setID&frame_ref==frame;
    nc_x_vec = nc_x_ref(frame_set_filter);
    nc_y_vec = nc_y_ref(frame_set_filter); 
    
    % indexing vectors    
    nc_sub_index_vec = sub_ind_ref(frame_set_filter); 
    nc_lin_index_vec = lin_ind_ref(frame_set_filter); 
    nc_master_vec = master_ind_ref(frame_set_filter);  
  
    % particle positions    
    pt_qc_vec = pt_qc_ref(frame_set_filter);
    spot_x_vec = spot_x_ref(frame_set_filter);
    spot_y_vec = spot_y_ref(frame_set_filter);        
    spot_z_vec = spot_z_ref(frame_set_filter);
    spot_x_vec3D = spot_x_ref3D(frame_set_filter);
    spot_y_vec3D = spot_y_ref3D(frame_set_filter);        
    spot_z_vec3D = spot_z_ref3D(frame_set_filter);
    particle_id_vec = pt_ref(frame_set_filter);
    src = set_key(set_key.setID==setID,:).prefix{1};

    % load stacks    
    mcp_stack = load_stacks(rawPath, src, frame, mcp_channel);
    protein_stack = load_stacks(rawPath, src, frame, proteinChannel);  
    
    [yDim, xDim, zDim] = size(protein_stack);    
    [x_ref,y_ref,z_ref] = meshgrid(1:xDim,1:yDim,1:zDim);
    % generate lookup table of inter-nucleus distances
    x_dist_mat = repmat(nc_x_vec,numel(nc_x_vec),1)-repmat(nc_x_vec',1,numel(nc_x_vec));
    y_dist_mat = repmat(nc_y_vec,numel(nc_y_vec),1)-repmat(nc_y_vec',1,numel(nc_y_vec));
    r_dist_mat = sqrt(double(x_dist_mat).^2 + double(y_dist_mat).^2);            
    
    % initialize arrays to store relevant info 
    for j = 1:numel(new_vec_fields)
        eval([new_vec_fields{j} ' = NaN(size(spot_x_vec));']);
    end
    for j = 1:numel(new_snip_fields)
        eval([new_snip_fields{j} ' = NaN(2*pt_snippet_size+1,2*pt_snippet_size+1,numel(spot_x_vec));']);
    end    
    
    % iterate through spots
    qc_mat = struct;
    for j = 1:numel(nc_x_vec)            
        % get location info
        x_nucleus = round(nc_x_vec(j));
        y_nucleus = round(nc_y_vec(j));   
        x_spot = round(spot_x_vec(j));
        y_spot = round(spot_y_vec(j));       
        z_spot = round(spot_z_vec(j))-1; % NL: why?
        % get position info from 3D fit
        x_spot3D = spot_x_vec3D(j);
        y_spot3D = spot_y_vec3D(j);
        z_spot3D = min([zMax spot_z_vec3D(j)-1]);        
        
        
        if isnan(x_spot) || ~pt_qc_vec(j)
            continue
        end
        % extract mask 
        spot_nc_mask = nc_ref_frame == nc_master_vec(j);                       
        %%%%%%%%%%%%%%%%%%%%%% Sample protein levels %%%%%%%%%%%%%%%%%%%%%%
        % get frames
        protein_frame = protein_stack(:,:,z_spot);
        mcp_frame = mcp_stack(:,:,z_spot);
        int_id = spot_roi_frame(y_spot,x_spot);        
        
        % regardless of qc issues filtered for later on, take pt samples
        % in vicinity of spot       
        spot_protein_vec(j) = nanmean(protein_frame(int_id==spot_roi_frame));
        spot_mcp_vec(j) = nanmean(mcp_frame(int_id==spot_roi_frame));      
        
        % volume protein sampling 
        spot_protein_vec_3d(j) = sample_protein_3D(x_spot3D,y_spot3D,z_spot3D,x_ref,y_ref,z_ref,xy_sigma,z_sigma,protein_stack);
        % make sure size is reasonable and that spot is inside nucleus
        if sum(spot_nc_mask(:)) < min_area || sum(spot_nc_mask(:)) > max_area || ~spot_nc_mask(y_spot,x_spot) %|| int_it==0  
            edge_qc_flag_vec(j) = -1;            
            serial_qc_flag_vec(j) = -1;
            continue
        end 
                      
        % sample snippets
        spot_protein_snips(:,:,j) = sample_snip(x_spot,y_spot,pt_snippet_size,protein_frame,spot_nc_mask);
        spot_mcp_snips(:,:,j) = sample_snip(x_spot,y_spot,pt_snippet_size,mcp_frame,spot_nc_mask);                 

        % Take average across all pixels within 1.5um of nuclues center 
        dist_mat = bwdist(~spot_nc_mask);        
        mf_samp_mask = dist_mat*PixelSize >= mf_samp_rad;              
        mf_null_protein_vec(j) = nanmean(protein_frame(mf_samp_mask));% / voxel_size;           
        % Edge sampling 
        spot_edge_dist = nc_dist_frame(y_spot,x_spot);        
        nc_edge_dist_vec = nc_dist_frame(spot_nc_mask);
        spot_edge_dist_vec(j) = spot_edge_dist;        
        spot_sep_vec = spot_dist_frame(spot_nc_mask);
        nc_indices = find(spot_nc_mask);
        
        % Now find control "spot" that is same distance from nucleus edge
        % as true spot
        [edge_null_x_vec(j), edge_null_y_vec(j), edge_null_nc_vec(j), edge_qc_flag_vec(j),~]...
            = find_control_sample(nc_edge_dist_vec, x_ref, y_ref, spot_sep_vec, spot_edge_dist,...
                 j, minSampleSep, spot_nc_mask,0);  
        % if initial attempt failed, try nearest neighbor nucleus
        null_mask = spot_nc_mask;
        if edge_qc_flag_vec(j) == 0
            % Find nearest neighbor nucleus
            r_vec = r_dist_mat(:,j);
            r_vec(j) = Inf;
            [~, mi] = min(r_vec);
            
            % get nn nucleus mask   
            x_spot_nn = round(spot_x_vec(mi));
            y_spot_nn = round(spot_y_vec(mi)); 
            
            nn_nc_mask = nc_ref_frame == nc_master_vec(mi);
            
            null_mask = nn_nc_mask; % reassign null mask
            if ~isnan(x_spot_nn)
                nan_flag = isnan(nn_nc_mask(y_spot_nn,x_spot_nn));
            end
            % make sure size is reasonable 
            if sum(nn_nc_mask(:)) >= min_area && sum(nn_nc_mask(:)) <= max_area                    
                nn_edge_dist_vec = nc_dist_frame(nn_nc_mask);
                nn_sep_vec = spot_dist_frame(nn_nc_mask);
                [edge_null_x_vec(j), edge_null_y_vec(j), edge_null_nc_vec(j), edge_qc_flag_vec(j),~]...
                    = find_control_sample(nn_edge_dist_vec, x_ref, y_ref, nn_sep_vec, spot_edge_dist,...
                         mi, minSampleSep, null_mask,1);
                if  edge_qc_flag_vec(j) == 1
                     edge_qc_flag_vec(j) =  2; % to flag cases when nn was used
                end
            end
        end 
        % Draw control samples (as appropriate)    
        if edge_qc_flag_vec(j) > 0  
            xc = edge_null_x_vec(j);
            yc = edge_null_y_vec(j);     

            null_dist_frame = zeros(size(protein_frame));
            null_dist_frame(yc,xc) = 1;
            null_dist_frame = bwdist(null_dist_frame);
            edge_null_protein_vec(j) = nanmean(protein_frame(nc_ref_frame>0&null_dist_frame<roi_rad_spot_pix));% / voxel_size;            
            edge_null_mcp_vec(j) = nanmean(mcp_frame(nc_ref_frame>0&null_dist_frame<roi_rad_spot_pix));% / voxel_size;            
            % take 3D protein sample
            edge_null_protein_vec_3d(j) = sample_protein_3D(xc,yc,z_spot3D,x_ref,y_ref,z_ref,xy_sigma,z_sigma,protein_stack);
            % draw snips    
            edge_null_protein_snips(:,:,j) = sample_snip(xc,yc,pt_snippet_size,protein_frame,nc_ref_frame>0);
            edge_null_mcp_snips(:,:,j) = sample_snip(xc,yc,pt_snippet_size,mcp_frame,nc_ref_frame>0);                                        
        end                  
        
        % Draw serialized control
        nc_index = nc_lin_index_vec(j);
        nc_sub_index = nc_sub_index_vec(j);   
        frame_vec_temp = nucleus_struct(nc_index).frames; 
        serial_null_x = nucleus_struct(nc_index).serial_null_x_vec;
        serial_null_y = nucleus_struct(nc_index).serial_null_y_vec;          
        % if this is the first sample for this spot, just find random
        % control snip. This will "seed" subsequent samples
        if all(isnan(serial_null_x))            
            % Now take a random sample                           
            sample_index_vec = 1:numel(spot_sep_vec);
            % filter for regions far enough away from locus
            cr_filter = spot_sep_vec >= minSampleSep & nc_edge_dist_vec >= minEdgeSep;
            sample_index_vec = sample_index_vec(cr_filter);
            % if candidate found, then proceed. Else look to neighboring nuclei
            if ~isempty(sample_index_vec)
                new_index = randsample(sample_index_vec,1);
                x_pos_vec = x_ref(spot_nc_mask);
                y_pos_vec = y_ref(spot_nc_mask);
                xc = x_pos_vec(new_index);
                yc = y_pos_vec(new_index);
                ec = nc_edge_dist_vec(new_index);                           
            else
                error('Unable to draw random sample. Check "PixelSize" and "min_sample_sep" variables')
            end  
        % otherwise, draw snip based on previous location
        else
            prev_frame = find(~isnan(serial_null_x),1,'last');
            n_frames = frame - frame_vec_temp(prev_frame); % used to adjust jump weights
            old_x = double(serial_null_x(prev_frame));
            old_y = double(serial_null_y(prev_frame));            
            % possible locations
            x_pos_vec = x_ref(spot_nc_mask);
            y_pos_vec = y_ref(spot_nc_mask);
            drControl = double(sqrt((old_x-x_pos_vec).^2+(old_y-y_pos_vec).^2));   
%             edge_dev_vec = spot_edge_dist-nc_edge_dist_vec;
            % calculate weights
            wt_vec = exp(-.5*((drControl/double((n_frames*driftTol/PixelSize))).^2));%+((edge_dev_vec)/roi_spot).^2));
            % anything too close to locus or with an edge distance too different from locus is excluded
            wt_vec(spot_sep_vec<minSampleSep|nc_edge_dist_vec < minEdgeSep) = 0;
            % draw sample
            xc = NaN;
            yc = NaN;
            ec = NaN; 
            if any(wt_vec>0)
                new_index = randsample(1:numel(x_pos_vec),1,true,wt_vec);
                xc = x_pos_vec(new_index);
                yc = y_pos_vec(new_index);
                ec = nc_edge_dist_vec(new_index);       
            else
                new_index = randsample(1:numel(x_pos_vec),1,true);
                xc = x_pos_vec(new_index);
                yc = y_pos_vec(new_index);
                ec = nc_edge_dist_vec(new_index);  
%             else
%                 error('This should not happen')
            end
        end       
        % draw samples
        null_dist_frame = zeros(size(protein_frame));
        if ~isnan(xc)
            null_dist_frame(yc,xc) = 1;
        end
        null_dist_frame = bwdist(null_dist_frame);
        % samples below default to NaN if no sample taken
        % sample protein
        serial_qc_flag_vec(j) = ~isnan(xc);
        serial_null_protein_vec(j) = nanmean(protein_frame(nc_ref_frame>0&null_dist_frame<roi_rad_spot_pix));% / voxel_size;        
        % record
        serial_null_x_vec(j) = xc;
        serial_null_y_vec(j) = yc;
        serial_null_edge_dist_vec(j) = ec;
        % 3D version                
        serial_null_protein_vec_3d(j) = sample_protein_3D(xc,yc,z_spot3D,x_ref,y_ref,z_ref,xy_sigma,z_sigma,protein_stack);% / sum(vol_denominator(:)) / voxel_size;                 
        % check for presence of sister spot
        x_spot_sister = NaN;
        y_spot_sister = NaN;
        if sum(nc_master_vec(j)==nc_master_vec) == 2
            indices = find(nc_master_vec(j)==nc_master_vec);
            si = indices(indices~=j);
            x_spot_sister = spot_x_vec(si);
            y_spot_sister = spot_y_vec(si);
        end            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        % save qc data                 
        qc_mat(numel(nc_x_vec)-j+1).setID = setID; 
        qc_mat(numel(nc_x_vec)-j+1).frame = frame;
        qc_mat(numel(nc_x_vec)-j+1).nc_index = nc_lin_index_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).nc_sub_index = nc_sub_index_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).qc_flag = edge_qc_flag_vec(j);            
        qc_mat(numel(nc_x_vec)-j+1).xp = x_spot;
        qc_mat(numel(nc_x_vec)-j+1).yp = y_spot;  
        qc_mat(numel(nc_x_vec)-j+1).xp_sister = x_spot_sister;
        qc_mat(numel(nc_x_vec)-j+1).yp_sister = y_spot_sister;  
        qc_mat(numel(nc_x_vec)-j+1).xc_edge = edge_null_x_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).yc_edge = edge_null_y_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).xc_serial = serial_null_x_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).yc_serial = serial_null_y_vec(j);       
        qc_mat(numel(nc_x_vec)-j+1).ParticleID = particle_id_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).serial_qc_flag = serial_qc_flag_vec(j);
        qc_mat(numel(nc_x_vec)-j+1).edge_qc_flag = edge_qc_flag_vec(j);        
        sz = nb_size;
        edge_dist_mat = nc_dist_frame;
        edge_dist_mat(~spot_nc_mask&~null_mask) = 0;
        if edge_qc_flag_vec(j) == 2         
            sz = max([nb_size,abs(x_nucleus - edge_null_x_vec(j)),abs(y_nucleus - edge_null_y_vec(j))...
                abs(x_nucleus - serial_null_x_vec(j)),abs(y_nucleus - serial_null_y_vec(j))]);            
        end
        y_range = max(1,y_nucleus-sz):min(yDim,y_nucleus+sz);
        x_range = max(1,x_nucleus-sz):min(xDim,x_nucleus+sz);
        qc_mat(numel(nc_x_vec)-j+1).x_origin = x_range(1);
        qc_mat(numel(nc_x_vec)-j+1).y_origin = y_range(1);
        qc_mat(numel(nc_x_vec)-j+1).mcp_snip = mcp_frame(y_range,x_range);
        qc_mat(numel(nc_x_vec)-j+1).protein_snip = protein_frame(y_range,x_range);
        qc_mat(numel(nc_x_vec)-j+1).edge_dist_snip = edge_dist_mat(y_range,x_range);             
    end 
    qc_structure(i).qc_mat = fliplr(qc_mat);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save snip data   
    % initialize struct to store snip data
    snip_data = struct;    
    % map data back to nucleus_struct    
    for j = 1:numel(nc_master_vec)
        nc_index = nc_lin_index_vec(j);
        nc_sub_index = nc_sub_index_vec(j);
        frame = nucleus_struct(nc_index).frames(nc_sub_index);
        for k = 1:numel(new_vec_fields)
            vec = eval(new_vec_fields{k});
            nucleus_struct(nc_index).(new_vec_fields{k})(nc_sub_index) = vec(j);
        end                             
        % store snips
        for k = 1:numel(new_snip_fields)
            snip = eval([new_snip_fields{k} '(:,:,j)']);                  
            snip_data.(new_snip_fields{k})(:,:,j) = snip;
        end  
    end
    
    % store key ID variables
    snip_data.frame = frame;
    snip_data.setID = setID;
    snip_data.particle_id_vec = particle_id_vec;
    snip_data.spot_edge_dist_vec = spot_edge_dist_vec;
    % indexing vectors    
    snip_data.nc_sub_index_vec = nc_sub_index_vec; 
    snip_data.nc_lin_index_vec = nc_lin_index_vec; 
    snip_data.nc_master_vec = nc_master_vec;    
    % specify name
     % read snip file    
    snip_name = ['snip_data_F' sprintf('%03d',frame) '_S' sprintf('%02d',setID)]; 
    if write_snip_flag            
        blank = struct;
        save([DataPath 'snip_data.mat'],'blank','-v7.3')    
        write_snip_flag = true;
    end
    snip_file = matfile([DataPath 'snip_data.mat'],'Writable',true);    
    snip_file.(snip_name)= snip_data;        
    clear snip_file;    
    % report time
    t = round(toc);
    disp([num2str(i) ' of ' num2str(size(set_frame_array,1)) ' frames completed (' num2str(t) ' sec)'])         
end
disp('saving qc frames...')
% save qc data
tic
particle_index = unique([nucleus_struct.ParticleID]);
particle_index = particle_index(~isnan(particle_index));
qc_particles = randsample(particle_index,min([100,numel(particle_index)]),false);
particle_index_full = [];
particle_frames_full = [];
for i = 1:numel(qc_structure)
    qc_mat = qc_structure(i).qc_mat;
    for  j = 1:numel(qc_mat)
        qc_spot = qc_mat(j);
        if ~isfield(qc_spot,'ParticleID')
            continue
        end
        ParticleID = qc_spot.ParticleID;
        if isempty(ParticleID) || ~ismember(ParticleID,qc_particles)
            continue
        end        
        frame = qc_spot.frame;      
        particle_index_full = [particle_index_full ParticleID];
        particle_frames_full = [particle_frames_full frame];        
        save_name = [snipPath 'pt' num2str(1e4*ParticleID) '_frame' sprintf('%03d',frame) '.mat'];
        save(save_name,'qc_spot');
    end
end
[particle_index_full, si] = sort(particle_index_full);
particle_frames_full = particle_frames_full(si);

qc_ref_struct.particle_frames_full = particle_frames_full;
qc_ref_struct.particle_index_full = particle_index_full;
toc
% save updated nucleus structure
disp('saving nucleus structure...')
nucleus_struct_protein = nucleus_struct;
save([DataPath 'qc_ref_struct.mat'],'qc_ref_struct')
save([DataPath 'nucleus_struct_protein.mat'],'nucleus_struct_protein','-v7.3') 