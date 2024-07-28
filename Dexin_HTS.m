addpath(genpath('F:\科研\博士后\胶原与预后\For Dexin\'));

clear all
close all
clc

filearray = dir('F:\刘文居\1-37-LSM格式\食管胃结合部腺癌--13-37\新建文件夹\*.lsm');
fea = zeros(length(filearray),142);
Fname = cell(length(filearray),1);
%load meta.mat;

for K=1:1:length(filearray)
    %if isnan(fea(K,7)) == 1

    fileName = ['F:\刘文居\1-37-LSM格式\食管胃结合部腺癌--13-37\新建文件夹\' filearray(K).name];
    fileName1 = ['F:\刘文居\1-37-LSM格式\食管胃结合部腺癌--13-37\新建文件夹\' filearray(K).name(1:end-3) 'tif'];
    disp(filearray(K).name(1:end-4));
    Fname{K,1} = filearray(K).name(1:end-4);
    %fileName = 'G:\early cancer 20171015\0801087-03.lsm';

    [Data, LSMinfo] = lsmread(fileName);
    width = LSMinfo.dimX;
    height = LSMinfo.dimY;
    channel1 = zeros(height,width);
    %channel2 = zeros(height,width);
    for i=1:1:height
        for j=1:1:width
            channel1(i,j) = Data(1,1,1,i,j);
            %channel2(i,j) = Data(1,2,1,i,j);
        end
    end
    if max(max(channel1)) > 255
        channel1 = uint8(channel1./2^4);
        channel2 = medfilt2(channel1,[3 3]);
        %collagen = otsu(channel1,2);
        %collagen_mask = collagen == 2;
        collagen_mask = channel2 > 30; %first batch 5
        %if sum(sum(collagen_mask))/(height*width) < 0.05
        %    collagen_mask = channel2 > 3;
        %end
    else
        channel2 = medfilt2(channel1,[3 3]);
        collagen_mask = channel2 > 80;
        %if sum(sum(collagen_mask))/(height*width) < 0.05
        %    collagen_mask = channel2 > 3;
        %end
    end
    collagen_mask = bwareaopen(collagen_mask,5);
    collagen_mask = imresize(collagen_mask,0.5);
    height = size(collagen_mask,1);
    width = size(collagen_mask,2);
    imwrite(collagen_mask,fileName1);
    if (sum(sum(collagen_mask))/(height*width)) < 0.01
        fea(K,1) = sum(sum(collagen_mask))/(height*width);
    else
    %morphology features
        p.Nimages   = 1;  
        p.yred      = 1:height;
        p.xred      = 1:width; 
        p = param_example(p); 
        im3 = zeros(1,height,width);
        im3(1,:,:) = collagen_mask;
        seg = zeros(1,height,width);
        seg(1,:,:) = logical(collagen_mask);
        data = fire(p,im3,2,seg);
        M = network_stat(data.Xa,data.Fa,data.Va,data.Ra);

        %M.area = sum(sum(collagen_mask))/sum(sum(mask));
        MMP = LSMinfo.voxSizeX*1000*1000/4;
        %fea(K,1) = M.area;
        fea(K,1) = sum(sum(collagen_mask))/(height*width);
        fea(K,2) = M.fiber_num/((height*width)*MMP*MMP);
        fea(K,3) = mean(M.L.*MMP);
        fea(K,4) = mean(M.RF.*MMP);
        fea(K,5) = mean(M.fstr);
        fea(K,6) = M.xlinkdens;
        fea(K,7) = mean(M.xlinkspace);
        [ori,~,~,~,F,~]=ftmethod2(channel1);
        fea(K,8) = min(F(1,1),F(2,2))/max(F(1,1),F(2,2));

        intensity features (first order histgram)
        Stats = chip_histogram_features( channel1,'NumLevels',16,'G',[] );
        fea(K,9:14) = Stats;

        %GLCM features
        text_fea = [];
        offsets = [0 1;-1 1;-1 0;-1 -1];
        for i=1:1:5
            offset = offsets.*i;
            for j=1:1:4
                glcm = graycomatrix(channel1,'Offset',offset(j,:));
               stats = graycoprops(glcm);
                text_fea = [text_fea stats.Contrast,stats.Correlation,stats.Energy,stats.Homogeneity];
            end
        end
        fea(K,15:94) = text_fea;

        %Gabor features
        gaborArray = gaborFilterBank(4,6,32,32);  % Generates the Gabor filter bank
        fea(K,95:142) = gaborFeatures_modified(channel1,gaborArray); 
        end
    end
    save('meta.mat','fea','Fname');
end
%save('meta.mat','fea','Fname');

