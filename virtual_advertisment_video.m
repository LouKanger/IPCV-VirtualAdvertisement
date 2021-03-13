%% L.W.J. Kanger - s1931318 - University of Twente
% Video from: https://www.youtube.com/watch?v=m-Vi_y0d-IE
% Match at stadium of AS Roma. Dimensions: 105 x 60 meters 
% (ref. https://www.football-italia.net/clubs/Roma/stadium)
%
% from 1.55 to 1.59/2.00 has clear view
%
% Additional rescourses
% https://github.com/cemunds/awesome-sports-camera-calibration

% Clear workspace and console. Close all figures
clear, clc 
close all

% Set debug to true to see the resulting image after each step
debug = false;

% Suppress specific warning
warning('off', 'Images:initSize:adjustingMag');

%% Load the video
filename = 'soccer_video_AJAX_1.mp4';
data_dir_name = 'data';
video_reader = VideoReader(data_dir_name+"/"+filename);

im_rgb_raw = readFrame(video_reader);
im_rgb = im2double(im_rgb_raw);
im = rgb2gray(im_rgb);

figure('units','normalized','position',[0.1 0.1 0.7 0.7]);
imshow(im_rgb,[],'InitialMagnification',200);

%% Create 2D exact field lines and show all intersection points
% International standards for the soccer field dimensions
width = 105; height = 68;
[field_lines, field_points] = generate_field_lines(width, height);
% first 7 field lines are vertical, rest is horizontal. num 4,6 lines left side 
field_lines_vert = field_lines(1:7,:);
field_lines_hor = field_lines(8:end,:);

% Create virtual advertisment sign in model coordinate system
w = 12;     % w meters long sign
h = 3;      % h meters high sign (flat on the ground)
d = 0.75;   % d meters from back field line
point1 = [-d-h, height/4 - w; -d, height/4 - w; -d-h, height/4 - w; -d-h, height/4];
point2 = [-d-h, height/4; -d, height/4; -d, height/4 - w; -d, height/4];
point1(:,2) = point1(:,2) + 5;
point2(:,2) = point2(:,2) + 5;
adv_points_model = table2struct(table(point1, point2));

% Loop over the vertical lines, calc intersect with all horizontal lines
field_line_int_points = zeros(size(field_lines_vert,1) * size(field_lines_hor,1), 2);
k = 1;
for i = 1:size(field_lines_vert,1)
    for j = 1:size(field_lines_hor,1)
        % using advantage of homogeneous coordinate system to calc intersection
        pt = cross(field_lines_vert(i,:), field_lines_hor(j,:));
        pt = pt ./ pt(end);

        % store intersection point
        field_line_int_points(k,:) = [pt(1), pt(2)];
        k = k + 1;
    end
end

% Select only the intersections on the left side of the field
field_line_int_points_left = zeros(4 * 6, 2);
k = 1;

% Loop over the vertical lines, calc intersect with all horizontal lines
for i = 1:4
    for j = 1:6
        % using advantage of homogeneous coordinate system to calc intersection
        pt = cross(field_lines_vert(i,:), field_lines_hor(j,:));
        pt = pt ./ pt(end);

        % store intersection point
        field_line_int_points_left(k,:) = [pt(1), pt(2)];
        k = k + 1;
    end
end

%% Define all parameters of the field line detection algorithm
dlambda = 30;       % number of color points around dominant color peak
tl = 128;           % luminace threshold
td = 20;            % difference threshold
tau = 10;           % line width assumption (twice this value)
num_peaks = 12;      % number of Hough peaks
fill_gap = 50;      % maximum gap between two linesegments (which still counts as 1 line)
min_length = 200;   % minimum length of the found lines
dtheta = 3.5;       % if angle difference bewteen 2 lines less than dtheta, remove 1
drho = 50;          % if distance between 2 lines less than drho, remove 1
tmax = 12;          % line has angle less than tmax is labeled horizontal
sigma_r = 6;        % max distance between white pixel and line
num_iterations = 3; % number of line refinement iterations

%% Perform the white field line detection algorithm
% Construct a field mask and enhance the result
field_mask_raw = construct_field_mask(im_rgb_raw, dlambda);
field_mask = enhance_field_mask(field_mask_raw);

% Apply the field mask to the rgb image
im_field_masked = im_rgb_raw .* field_mask;
im_masked_gray = im2double(rgb2gray(im_field_masked));

% Detect the white field line pixels and apply the mask (for debuggin/visualisation)
white_field_lines = detect_white_pixels(im_field_masked, tl, td, tau);
im_white_masked = im_field_masked .* uint8(abs(double(white_field_lines) - 1));

% Perform a Hough transform on the binary lines/edges image to get field lines
lines = hough_line_detection(white_field_lines, num_peaks, fill_gap, min_length);

% Remove duplicate lines (line segments that represent the same field line)
lines = remove_duplicate_lines(lines, dtheta, drho);

% Divide lines into vertical and horizontal groups, sort them and convert to homogeneous coordinates
[lines_vert, lines_hor] = group_lines(lines, tmax);

% Detect and remove any outliers (goal post lines or arc lines)
angles_vert = [lines_vert.theta];
angles_hor = [lines_hor.theta];
lines_vert(isoutlier(angles_vert)) = [];
lines_hor(isoutlier(angles_hor)) = [];

% Convert to homogeneous coordinates
lines_vert_h = lines_to_homogeneous(lines_vert);
lines_hor_h = lines_to_homogeneous(lines_hor);

% Refine the line parameters
[lines_vert_h, lines_hor_h] = refine_line_parameters(white_field_lines, ...
    lines_vert_h, lines_hor_h, sigma_r, num_iterations);

% Calculate the intersections
intersections = zeros(size(lines_vert_h,1)*size(lines_hor_h,1), 2);
k = 1;
for i = 1:size(lines_vert_h,1)
    for j = 1:size(lines_hor_h,1)
        % calc intersection of two lines in homogeneous coordinates [cx,cy,c]
        pt = cross(lines_vert_h(i,:), lines_hor_h(j,:));
        
        % convert to Euclidian coordinates [x,y] and store point
        intersections(k,:) = [pt(1), pt(2)] ./ pt(3);
        k = k + 1;
    end
end

%% Manually select corresponding intersections
% create figure and show the detected field lines and intersections
figure('units','normalized','position',[0.1 0.1 0.7 0.7]);
suptitle('Write down the corresponding indices in the array')
subplot(2,1,1);
hold on
set(gca,'Ydir','reverse')
daspect([1,1,1])
scatter(intersections(:,1), intersections(:,2), 30, 'red', 'filled');
imshow(im_rgb_raw,[],'InitialMagnification',200);
xlim([0, size(im_rgb_raw,2)])
ylim([0, size(im_rgb_raw,1)])
hold on
plot_hlines(lines_vert_h, white_field_lines, 'blue');
plot_hlines(lines_hor_h, white_field_lines, 'blue');
scatter(intersections(:,1), intersections(:,2), 30, 'red', 'filled');
for i = 1:size(intersections, 1)
    pt = intersections(i,:);
    text(pt(1), pt(2)+30, num2str(i),'FontSize',12);
end

% show the model field lines and intersections
subplot(2,1,2);
hold on
set(gca,'Ydir','reverse')
axis off
xlim([-5, width/2 + 5])
daspect([1,1,1])
for k = 1:4     % vertical lines left side of the field
    xy = [field_points(k).point1; field_points(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
end
for k = 8:13    % horizontal lines left side of the field
    xy = [field_points(k).point1; field_points(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
end
scatter(field_line_int_points_left(:,1), field_line_int_points_left(:,2), 30, 'red', 'filled');
for i = 1:size(field_line_int_points_left, 1)
    pt = field_line_int_points_left(i,:);
    text(pt(1)+0.5, pt(2)+2, num2str(i),'FontSize',12);
end

% ask to write down the correct combination of intersections
point_idx = input('Write down the corresponding intersection points: [p1, q1; p2, q2; ...]\n');
% soccer_field_1: [1,2 ; 2,4 ; 3,5 ; 4,8 ; 5,10 ; 6,11 ; 7,14 ; 8,16; 9,17]
% soccer_field_3: [7,2 ; 8,4 ; 9,5 ; 10,8 ; 11,10 ; 12,11 ; 13,14 ; 14,16; 15,17]
% soccer_field_4: [1,2 ; 2,4 ; 3,5 ; 4,6 ; 5,8 ; 6,10 ; 7,11 ; 8,12]
% video_1: [11,2; 10,4; 9,5; 12,6; 15,8; 14,10; 13,11; 16,12]
% screen_shot_AJAX_1: [1,1; 2,2; 3,3; 4,4; 6,7; 7,8; 8,9; 9,10; 11,13; 12,14; 13,15; 14,16; 15,17]
% v2: [1,1; 2,2; 3,3; 4,4; 6,7; 7,8; 8,9; 12,14; 13,15; 14,16; 15,17]

% store the result (p = detected intersections, q = model intersections)
p_init = intersections(point_idx(:,1),:);
q_init = field_line_int_points_left(point_idx(:,2),:);

%% Determine the homography to warp advertisment sign onto the image
% Apply scale and translation to model for better result
scale = 10;
offset = (size(im,1) - scale * height) / 2;

% estimate the geometric transform using the hand-picked point pairs
Hm2w = estimateGeometricTransform(q_init * scale + offset, p_init, 'projective', ...
    'MaxNumTrials', 2000, 'Confidence', 99.9);

if debug == true
    % transform model field lines to the image coordinate frame
    % (transformPointsInverse is also available)
    points1 = zeros(length(field_points),2);
    points2 = zeros(length(field_points),2);
    for k = 1:length(points1)
        points1(k,:) = [field_points(k).point1] * scale + offset;
        points2(k,:) = [field_points(k).point2] * scale + offset;
    end
    world_points1 = transformPointsForward(Hm2w, points1);
    world_points2 = transformPointsForward(Hm2w, points2);
    field_points_world = table2struct(table(world_points1, world_points2, ...
        'VariableNames', {'point1', 'point2'}));

    % show original image with field lines projected on them
    figure('units','normalized','position',[0.1 0.1 0.7 0.7]);
    imshow(im_rgb_raw,[],'InitialMagnification',200);
    hold on
    for k = 1:length(field_points_world)
        xy = [field_points_world(k).point1; field_points_world(k).point2];
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');
    end
end

% Warp advertisment sign from model coordinates to image/world coordinates
points1 = zeros(length(adv_points_model),2);
points2 = zeros(length(adv_points_model),2);
for k = 1:length(points1)
    points1(k,:) = [adv_points_model(k).point1] * scale + offset;
    points2(k,:) = [adv_points_model(k).point2] * scale + offset;
end
world_points1 = transformPointsForward(Hm2w, points1);
world_points2 = transformPointsForward(Hm2w, points2);
adv_points_world = table2struct(table(world_points1, world_points2, ...
    'VariableNames', {'point1', 'point2'}));

% % show original image with field lines projected on them
% figure('units','normalized','position',[0.1 0.1 0.7 0.7]);
% imshow(im_rgb_raw,[],'InitialMagnification',200);
% hold on

X = [adv_points_world(1).point1(1), adv_points_world(2).point1(1), ...
    adv_points_world(2).point2(1), adv_points_world(1).point2(1)];
Y = [adv_points_world(1).point1(2), adv_points_world(2).point1(2), ...
    adv_points_world(2).point2(2), adv_points_world(1).point2(2)];

% patch(X,Y,'red');
% title('Frame 1');

%% Place the advertisment sign onto all frames of the video

% TODO:
% - As a final step, make the algorithm fully automatic. So remove the manual
%   selection of corresponding intersection points by finding the best fitting
%   homography (see the article that also does the white line pixel detection)

% Algorithm Steps:
% 0) Before algorithm starts: detect lines using the self-made algorithm and
%    manually select intersections in the image and the corresponding points in
%    the model. This determines the initial state and the points that will be
%    used in all the future frames.
% 1) Create white pixel mask of current frame
% 2) Refine the lines of the previous frame such that they match the new frame
% 3) Calcule the homography
% 4) Place advertisment sign
% 5) Repeat

% Load advertisment sign image
[ad_img, cmap] = imread(data_dir_name+"/"+"advertisment_image_1.png");
ad_img = im2uint8(ind2rgb(ad_img, cmap));
ad_img = imrotate(ad_img, 90);  % Rotate to place along back field line

% Setup for warping the advertisment sign
x = zeros(length(adv_points_model), 2);
y = zeros(length(adv_points_model), 2);
for k = 1:length(adv_points_model)
    pt1 = adv_points_model(k).point1;
    pt2 = adv_points_model(k).point2;
    x(k,1) = pt1(1) * scale + offset;
    x(k,2) = pt2(1) * scale + offset;
    y(k,1) = pt1(2) * scale + offset;
    y(k,2) = pt2(2) * scale + offset;
end

% Create the image reference object for the image in the scale model coordinates
xWorldLimits = [min(x, [], 'all'), max(x, [], 'all')];
yWorldLimits = [min(y, [], 'all'), max(y, [], 'all')];
ad_img_ref = imref2d(size(ad_img), xWorldLimits, yWorldLimits);

% Reset the video reader and get the number of frames (quirky Matlab stuff)
video_reader = VideoReader(data_dir_name+"/"+filename);
num_frames = video_reader.NumberOfFrames;
video_reader = VideoReader(data_dir_name+"/"+filename);

% Create video multi-frame matrix to store the frames
max_frames = 150; %30;

% Store the previous lines
lines_vert_h_old = lines_vert_h;
lines_hor_h_old = lines_hor_h;

% Array to store the corners of the transformed virtual advertisment sign
warped_ad_points = zeros(max_frames, 2, 4); % (num_frames, xy, 4 corners)

% Frame counter
frame = 1;

disp("Transforming advertisment sign onto the frames.")
tic;
while hasFrame(video_reader)
    disp("Frame: " + num2str(frame))
    
    % Grab next frame of the video reader
    im_rgb_raw = readFrame(video_reader); 
    im_rgb = im2double(im_rgb_raw);
    im = rgb2gray(im_rgb);
    
    % ------------------------------------------------------------------------ %
    % ------------[ Detect the white field lines in the new frame ]----------- %
    % ------------------------------------------------------------------------ %
    % Construct a field mask and enhance the result
    field_mask_raw = construct_field_mask(im_rgb_raw, dlambda);
    field_mask = enhance_field_mask(field_mask_raw);

    % Apply the field mask to the rgb image
    im_field_masked = im_rgb_raw .* field_mask;

    % Detect the white field line pixels
    white_field_lines = detect_white_pixels(im_field_masked, tl, td, tau);
    
    % ------------------------------------------------------------------------ %
    % ---------------[ Find the field lines and the key POIs ]---------------- %
    % ------------------------------------------------------------------------ %
    % Find the new lines by refining the old lines using the new frame
    [lines_vert_h, lines_hor_h] = refine_line_parameters(white_field_lines, ...
        lines_vert_h_old, lines_hor_h_old, sigma_r, num_iterations);
    
    % Update the intersections and lines
    lines_vert_h_old = lines_vert_h;
    lines_hor_h_old = lines_hor_h;

    % Calculate the intersections of the new lines
    intersections = zeros(size(lines_vert_h,1)*size(lines_hor_h,1), 2);
    k = 1;
    for i = 1:size(lines_vert_h,1)
        for j = 1:size(lines_hor_h,1)
            % calc intersection of two lines in homogeneous coordinates [cx,cy,c]
            pt = cross(lines_vert_h(i,:), lines_hor_h(j,:));

            % convert to Euclidian coordinates [x,y] and store point
            intersections(k,:) = [pt(1), pt(2)] ./ pt(3);
            k = k + 1;
        end
    end

    % ------------------------------------------------------------------------ %
    % ----------[ Calculate Homography and transform advertisment ]----------- %
    % ------------------------------------------------------------------------ %
    % Get the matching points using the new intersections
    p = intersections(point_idx(:,1),:);
    q = field_line_int_points_left(point_idx(:,2),:);  % constant
    
    % Estimate the geometric transform using the point pairs
    Hm2w = estimateGeometricTransform(q * scale + offset, p, 'projective', ...
        'MaxNumTrials', 2000, 'Confidence', 99.9);

    % Warp advertisment sign from model coordinates to image/world coordinates
    points1 = zeros(length(adv_points_model),2);
    points2 = zeros(length(adv_points_model),2);
    for k = 1:length(points1)
        points1(k,:) = [adv_points_model(k).point1] * scale + offset;
        points2(k,:) = [adv_points_model(k).point2] * scale + offset;
    end
    world_points1 = transformPointsForward(Hm2w, points1);
    world_points2 = transformPointsForward(Hm2w, points2);
    adv_points_world = table2struct(table(world_points1, world_points2, ...
        'VariableNames', {'point1', 'point2'}));

    % Store the transformed advertisment coordinates
    X = [adv_points_world(1).point1(1), adv_points_world(2).point1(1), ...
        adv_points_world(2).point2(1), adv_points_world(1).point2(1)];
    Y = [adv_points_world(1).point1(2), adv_points_world(2).point1(2), ...
        adv_points_world(2).point2(2), adv_points_world(1).point2(2)];
    
    warped_ad_points(frame, 1, :) = X(:);
    warped_ad_points(frame, 2, :) = Y(:);
   
    % Update the frame counter and check if max frame number is reached
    frame = frame + 1;
    if frame > max_frames
        break
    end
end

calc_time = toc;
disp("Completed calculations for each frame in " + calc_time + " seconds.")

%% Create the video
% Smooth the transformed virtual advertisment points
frame_span = 15;
ad_pos = zeros(size(warped_ad_points, 1), 4, 2);
for i = 1:size(ad_pos,2)
    ad_pos(:,i,1) = smooth(warped_ad_points(:,1,i), frame_span, 'rloess');
    ad_pos(:,i,2) = smooth(warped_ad_points(:,2,i), frame_span, 'rloess');
end

% Define the new image reference object for the advertisment
qq = [0,                0
      size(ad_img,1),   0;
      size(ad_img,1),   size(ad_img,2);
      0,                size(ad_img,2)];

qq_xWorldLimits = [min(qq(:,1)), max(qq(:,1))];
qq_yWorldLimits = [min(qq(:,2)), max(qq(:,2))];
ad_img_ref2 = imref2d(size(ad_img), qq_xWorldLimits, qq_yWorldLimits);

% Reset the video reader
video_reader = VideoReader(data_dir_name+"/"+filename);

% Create video object with 15 FPS
video_name = 'ajaxVidVirtualAd_v5.avi';
writerObj = VideoWriter(video_name);
writerObj.FrameRate = 15;

disp("Creating video ...")

% Write the frames to the video object to create the video
open(writerObj);
fig = figure('units','normalized','position',[0.1 0.1 0.7 0.7]);
imshow(im_rgb_raw,[],'InitialMagnification',200);
hold on

tic;
for frame = 1:size(ad_pos, 1)
    frame_img = readFrame(video_reader);
    
    % Find the geometric transform of the ad image to the virtual position
    Tform = estimateGeometricTransform(qq, reshape(ad_pos(frame,:,:), 4,2), ...
        'projective', 'MaxNumTrials', 2000, 'Confidence', 99.9);
    
    % Warp the advertisment image onto the field
    warped_ad_img = imwarp(ad_img, ad_img_ref2, Tform, ...
        'OutputView', imref2d(size(im_rgb_raw)));

    % Place only the advertisment sign on the green field pixels of the new frame
    field_mask = construct_field_mask(frame_img, round(dlambda*1.25));
    field_mask = field_mask .* ones([size(field_mask),3]);
    masked_warped_ad = im2uint8(im2double(warped_ad_img) .* field_mask);

    idx = (sum(masked_warped_ad, 3) ~= 0);
    idx = logical(idx .* ones([size(idx),3]));
    stitched = frame_img;
    stitched(idx) = masked_warped_ad(idx);
    
    imshow(stitched,[],'InitialMagnification',200);
%     patch(ad_pos(frame,:,1), ad_pos(frame,:,2),'red');
    title("Frame " + num2str(frame));
    
    % Store the figure frame
    writeVideo(writerObj, getframe(fig));
    
    % Clear figure data to reduce memory usage
    clf(fig);
end
close(writerObj);

vid_create_time = toc;
disp("Finished video with name: '"+video_name+"'")
disp("Calculation time is " + vid_create_time + " seconds.")

%% Plot the corner positions of the virtual advertisment sign (data visualisation)
figure;
ax1 = subplot(1,2,1);
hold on
for i = 1:4
    p1 = warped_ad_points(:, :, i);
    h = plot(p1(:,1), '.', 'MarkerSize', 5);
    plot(smooth(p1(:,1), 0.1, 'rloess'), 'Color', h.Color, 'LineWidth', 1.5);
end
hold off
title("x pos");
legend("Raw 1", "Filtered 1", "Raw 2", "Filtered 2", "Raw 3", "Filtered 3", ...
    "Raw 4", "Filtered 4", 'Location', 'southeast', 'Position',[0.35 0.12 0.17 0.28]);
xlabel("Frame number");
ylabel("Pixel number");

ax2 = subplot(1,2,2);
hold on
for i = 1:4
    p1 = warped_ad_points(:, :, i);
    h = plot(p1(:,2), '.', 'MarkerSize', 5);
    plot(smooth(p1(:,2), 0.1, 'rloess'), 'Color', h.Color, 'LineWidth', 1.5);
end
hold off
title("y pos");
% legend("Raw 1", "Filtered 1", "Raw 2", "Filtered 2", "Raw 3", "Filtered 3", ...
%     "Raw 4", "Filtered 4", 'Location', 'northwest');
xlabel("Frame number");
ylabel("Pixel number");



