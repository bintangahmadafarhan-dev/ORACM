clear; close all; clc;
list = dir('images\'); list = list(4:end,1);
show = 1; result = zeros(0);
ImageIndex = 1; % SELECT Image from (1-15)
if ~strcmp(list(ImageIndex,1).name,{'.','..','desktop.ini'})        
    img = imread(['images\',list(ImageIndex,1).name]);
    if size(img,3)>1,img = rgb2gray(img);end
    img = double(img);
    [time1,itr1] = ORACM(img,show);    
    [time2,itr2] = ACMwithSBGFRLS(img,show);   
    [time3,itr3] = CHENVESE(img,show);
end