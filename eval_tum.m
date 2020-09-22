clear all;
close all;

path = '/Users/yyqing/Soft/Project/gitOrb3/orb_slam3/';
dataset_path = '/Users/yyqing/Soft/DataSet/tum/dataset-outdoors5_512_16/mav0/mocap0/data.csv';

trajectory = 'kf_dataset-corridor1_512_monoi.txt';

file_name = [path,trajectory];

fprintf(file_name);

[timestamp, x, y, z, q0, q1, q2, q3] = textread(file_name);

figure;
plot(x, y, 'r.');
title('est trajectory');

% fid = fopen(dataset_path);
% Title = textscan(fid, '%s %s %s %s %s %s %s %s',1,'delimiter', ',');
% Data = textscan(fid, '%f %f %f %f %f %f %f %f','delimiter', ',');
% grouth_x = cell2mat(Data(2));
% grouth_y = cell2mat(Data(3));
% fclose(fid);
% 
% figure;
% plot(grouth_x, grouth_y, 'g.');
% title('groud truth trajectory');