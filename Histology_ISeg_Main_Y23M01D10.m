%
% Code to manually segment histology section images and quantify stain within atherosclerotic plaques.
%  
%
% R. Pedrigi (Dec, 2018)


tic;

%% Inputs

% Input the directory of the root folder that contains all mouse histology files (seperated into folders, one per mouse).
rootdir = ['F:\FlowRawData\Histology\'];               % Example: C:\myfiles\Histology\

% Input stain
stain = 'ORO';              % ORO for oil red O (lipid), PAN for CD68 (pan macrophages), and PSR for picrosirius red (collagen).

% PSR color to analyze
psrcol = 1;                 % PSR red and green were binarized separately and must be quantified separately; == 1 for red and == 2 for green

% Manual contouring
recontour = 1;              % == 1 for extracting stain from the new image using previously defined contours, == 0 for re-contouring.

% Save MAT file 
saveall = 0;                % == 1 to save a MAT file of quantified histology stains, == 0 to not save           
 
% Mouse number
mousenum = [1];


%% Main Code

for k = 1:size(mousenum,2)

    clc;
    clearvars -except 'mousenum' 'k' 'rootdir' 'stain' 'recontour' 'saveall' 'imagergb' 'psrcol';
    
    %%% Main code

    % Idenfity sub-directory name that contains histology image files for the specified mouse.
    if mousenum(k) < 10
        mouseid = strcat('M00',num2str(mousenum(k)));
    else
        mouseid = strcat('M0',num2str(mousenum(k)));
    end
    filedir1 = strcat(rootdir,mouseid,'\',stain,'\')
    
    % Identify specific color of PSR stain to analyze
    if strcmp(stain,'PSR') == 1
        if psrcol == 1
            datadir1 = strcat(filedir1,'RED\');
        else
            datadir1 = strcat(filedir1,'GREEN\');
        end
    else
        datadir1 = filedir1;
    end
    
    if exist(filedir1,'file') == 0
        disp('Folder does not exist');
        continue;
    else
        filedirimage1 = Histology_ISeg_Sub1_ArrayReshape_Y23M01D10(filedir1,datadir1,mouseid,stain,'M*');
    end

    switch stain

        case {'ORO','PSR'}      % Binary images to quantify with brightfield images for contouring
            filedirimage2 = Histology_ISeg_Sub1_ArrayReshape_Y23M01D10(filedir1,datadir1,mouseid,stain,'Binary*');
            imagenum = size(filedirimage2,1);
        case 'PAN'              % Fluorescence images to quantify with phase images for contouring
            filedirimage2 = Histology_ISeg_Sub1_ArrayReshape_Y23M01D10(filedir1,datadir1,mouseid,stain,'Phase*');
            imagenum = size(filedirimage2,1);
    end

    Econtour1 = cell(imagenum,1);
    Econtour2 = cell(imagenum,1);
    Econtour3 = cell(imagenum,1);

    pball = nan(40,2);
    stainlall = nan(40,4);
    stainpall = nan(40,4);
    stainwall = nan(40,5);
    
    % Loop through all histology images in sub-directory
    for i = 1:imagenum
        
        % Identify the directory and file name of each histological section.
        if strcmp(stain,'PSR') == 1
            strloc = strfind(filedirimage2{i,1},'Binary_');
            secfiledir1 = strcat(filedirimage2{i,1}(1:strloc-1),filedirimage2{i,1}(strloc+7:end-4));
            disp(filedirimage2{i,1}(strloc+7:end-4));
        else
            secfiledir1 = filedirimage1{i,:}(1:size(filedirimage1{i,:},2)-4);
            disp(filedirimage1{i,:});
        end

        listinglocal = dir(strcat(filedirimage1{i,:}));

        filename = listinglocal(1).name;

        %%% Load images
        switch stain

            case {'ORO','PSR'}
                E1 = imread(strcat(filedirimage1{i,:}));
                E2 = imread(strcat(filedirimage2{i,:}));

                Efinal = E2;

            case 'PAN'
                E1 = imread(strcat(filedirimage2{i,:}));
                E2 = imread(strcat(filedirimage1{i,:}));

            	Efinal = rgb2gray(E2);

        end
        imscale = 255;

        %%% Contour images
        if recontour == 0
            % Draw contours onto the phase image to isolate the vessel wall. You will be prompted to draw for the number of "imfreehand" commands below.
            disp('Draw three contours in the following order: (1) lumen, (2) intima (IEL), (3) EEL'); 
            prompt = 'Keep contour (==1), recontour (==0), or skip section (==2)?';

            itr = 1;
            userinput = 0;
            while userinput == 0

                if itr > 1
                    close all;
                end

                figure1 = figure('PaperType','A','Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1],'Name',filename);
                hold on;
                Ehandle = imshow(E1,'InitialMagnification','fit');
                disp(strcat(['Image',' ',num2str(i),' of',' ',num2str(imagenum)]));

                Econtour1{i,1} = imfreehand(gca,'Closed','True');                   %Draw lumen (yellow)
                setColor(Econtour1{i,1},[1 1 0]);

                Econtour2{i,1} = imfreehand(gca,'Closed','True');                   %Draw IEL (pink)
                setColor(Econtour2{i,1},[1 0 0.5]);

                Econtour3{i,1} = imfreehand(gca,'Closed','True');                   %Draw EEL (green)
                setColor(Econtour3{i,1},[0 1 0]);  

                userinput = input(prompt);    

                itr = itr+1;
            end

            % Obtain the xy-coordinates of the EEL and IEL contours
            Econtour1xy = Econtour1{i,1}.getPosition;
            Econtour2xy = Econtour2{i,1}.getPosition;
            Econtour3xy = Econtour3{i,1}.getPosition;

        else

         	load(strcat(secfiledir1,'.mat'),'Econtour1xy','Econtour2xy','Econtour3xy');

            if size(Econtour1xy,1) > 30
                userinput = 1;
            else
                userinput = 0;
            end

        end

        switch stain

            case {'ORO','PSR'}

                mask1 = poly2mask(Econtour1xy(:,1),Econtour1xy(:,2),size(Efinal,1),size(Efinal,2));
                mask2 = poly2mask(Econtour2xy(:,1),Econtour2xy(:,2),size(Efinal,1),size(Efinal,2));
                if size(Econtour3xy,1) < 20
                    mask3 = nan;
                else
                    mask3 = poly2mask(Econtour3xy(:,1),Econtour3xy(:,2),size(Efinal,1),size(Efinal,2));
                end

                Efinalinv = 255-Efinal;

                Emask0q = nan;                                  % Emask0 = endothelium

                Emask1 = Efinal.*0;                             % Emask1 = plaque.
                Emask1(mask2==1) = Efinalinv(mask2==1);         % Fill the inside of mask 2 (IEL) with the pixel values from image Efinal.
                Emask1(mask1==1) = 0;                           % Fill the inside of mask 1 (Lumen) with white (pixel value = 255).
                Emask1q = logical(Emask1);

                Emask2 = Efinal.*0;                             % Emask2 = vessel wall.
                Emask2(mask3==1) = Efinalinv(mask3==1);        	% Fill the inside of mask 3 (EEL) with the pixel values from image Efinal.
                Emask2(mask2==1) = 0;                           % Fill the inside of mask 2 (IEL) with white (pixel value = 255).
                Emask2q = logical(Emask2);

                Emask3q = nan;                                  % Lumen

            case 'PAN'

                mask1 = poly2mask(Econtour1xy(:,1),Econtour1xy(:,2),size(Efinal,1),size(Efinal,2));
                mask2 = poly2mask(Econtour2xy(:,1),Econtour2xy(:,2),size(Efinal,1),size(Efinal,2));
                if size(Econtour3xy,1) < 20
                    mask3 = nan;
                else
                    mask3 = poly2mask(Econtour3xy(:,1),Econtour3xy(:,2),size(Efinal,1),size(Efinal,2));
                end

                Emask1 = Efinal.*0;                                                              	% Emask1 = entire plaque
                Emask1(mask2==1) = Efinal(mask2==1);                                                    
                Emask1(mask1==1) = 0;                                                                   
                Emask1q = double(Emask1)./imscale;

                Emask2 = Efinal.*0;                                                              	% Emask2 = vessel wall.
                Emask2(mask3==1) = Efinal(mask3==1);                                              	% Fill the inside of mask 3 (EEL) with the pixel values from image Efinal.
                Emask2(mask2==1) = 0;                                                             	% Fill the inside of mask 2 (IEL) with white (pixel value = 255).
                Emask2q = double(Emask2)./imscale;

                Emask3q = nan;                                                                      % Emask3 = lumen

        end

        if userinput == 1

            pb = (sum(sum(mask2))-sum(sum(mask1)))/sum(sum(mask3));                                 % PB (plaque burden) = plaque area (IEL-lumen) to EEL area (including lumen)

            stainl = sum(sum(Emask3q))/sum(sum(mask1));                                          	% Lumen stain normalized by lumen area
            stainp = sum(sum(Emask1q))/(sum(sum(mask2))-sum(sum(mask1)));                           % Stain in the plaque of the vessel normalized by plaque area (or percent of plaque area stained because the stain is binarized)
            stainw = sum(sum(Emask2q))/(sum(sum(mask3))-sum(sum(mask2)));                           % Stain in the wall of the vessel normalized by wall area   

            stainlarea = sum(sum(Emask3q));                                                         % Stain area of the lumen
            stainparea = sum(sum(Emask1q));                                                         % Stain area in the plaque
            stainwarea = sum(sum(Emask2q));                                                         % Stain area in the wall

            larea = sum(sum(mask1));                                                                % Lumen area
            parea = (sum(sum(mask2))-sum(sum(mask1)));                                              % Plaque area
            warea = (sum(sum(mask3))-sum(sum(mask2)));                                              % Artery wall area
            sarea = sum(sum(mask3));                                                                % Area of the entire section (i.e., EEL area), including the lumen

        else

            pb = nan;                                

            stainl = nan;                                          	
            stainp = nan;                           
            stainw = nan;                           

            stainlarea = nan;                                                         
            stainparea = nan;                                                         
            stainwarea = nan;                                                       

            larea = nan;                                                                
            parea = nan;                                             
            warea = nan;                                              
            sarea = nan;                                                               

        end

        % Save variables to a MAT file specific to the vessel section
        if saveall == 1
            save(strcat(secfiledir1,'.mat'),'Econtour1','Econtour2','Econtour3','Econtour1xy','Econtour2xy','Econtour3xy',...
                        'mask1','mask2','mask3','Emask1','Emask1q','Emask2','Emask2q','stainl','stainlarea','larea',...
                        'stainp','stainparea','parea','stainw','stainwarea','warea','sarea','pb');
        end

        if contains(filename,strcat(filename(1:9),'_SI_')) == 1
            tmpnum = 0;
        elseif contains(filename,strcat(filename(1:9),'_SII_')) == 1
            tmpnum = 10;
        elseif contains(filename,strcat(filename(1:9),'_SIII_')) == 1
            tmpnum = 20;
        elseif contains(filename,strcat(filename(1:9),'_SIV_')) == 1
            tmpnum = 30;
        elseif contains(filename,strcat(filename(1:9),'_SV_')) == 1
            tmpnum = 40;
        elseif contains(filename,strcat(filename(1:9),'_SVI_')) == 1
            tmpnum = 50;
        end

        if isempty(pb) == 0
            pball(i,:) = [str2double(filename(end-5:end-4))+tmpnum pb]; 
        else 
            pball = [];
        end
        stainlall(i,:) = [str2double(filename(end-5:end-4))+tmpnum stainl stainlarea larea]; 

        stainpall(i,:) = [str2double(filename(end-5:end-4))+tmpnum stainp stainparea parea];

        stainwall(i,:) = [str2double(filename(end-5:end-4))+tmpnum stainw stainwarea warea sarea];

    end

    % Save global variables (over all sections) to a MAT file for the folder (mouse and stain).
    if saveall == 1
          save(strcat(datadir1,mouseid,'L_',stain,'.mat'),'pball','stainlall','stainpall','stainwall');
    end

end

toc;





