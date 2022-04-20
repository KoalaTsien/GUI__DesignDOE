%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright @ Koala
% 2020-10-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = GUI_DesignDOE()

close all

S.fh = figure('units','norm',...
    'position',[0.05 0.05 0.9 0.9],...
    'menubar','none',...
    'name','GS相位恢复算法 by Yanjun',...
    'numbertitle','off',...
    'resize','off');%,...
% 'WindowButtonDownFcn','disp(''gcf clicked！'')');

%% 波长
S.tx_wavelength = uicontrol('style','text',...
    'unit','norm',...
    'position',[0.01 0.9 0.15 0.03],...
    'string','波长/nm = ',...
    'fontsize',16);
S.ed_wavelength = uicontrol('style','edit',...
    'unit','norm',...
    'position',[0.15 0.9 0.05 0.03],...
    'string','600',...
    'fontsize',12);

%% 焦距
S.tx_f = uicontrol('style','text',...
    'unit','norm',...
    'position',[0.01 0.8 0.15 0.03],...
    'string','焦距/cm = ',...
    'fontsize',16);
S.ed_f = uicontrol('style','edit',...
    'unit','norm',...
    'position',[0.15 0.8 0.05 0.03],...
    'string','10',...
    'fontsize',12);

%% DOE像素数目
S.tx_NN = uicontrol('style','text',...
    'unit','norm',...
    'position',[0.01 0.7 0.15 0.03],...
    'string','DOE像素数目=',...
    'fontsize',16);
S.ed_NN = uicontrol('style','edit',...
    'unit','norm',...
    'position',[0.15 0.7 0.05 0.03],...
    'string','500',...
    'fontsize',12);

%% DOE像素边长
S.tx_pixeld = uicontrol('style','text',...
    'unit','norm',...
    'position',[0.05 0.51 0.15 0.03],...
    'string','DOE像素边长/um',...
    'fontsize',16);
S.ed_pixeld = uicontrol('style','edit',...
    'unit','norm',...
    'position',[0.05 0.47 0.1 0.03],...
    'string','5.00000',...
    'fontsize',12);
align([S.tx_pixeld S.ed_pixeld],'Center','Fixed',7);




%% 成像面总边长
S.tx_fL = uicontrol('style','text',...
    'unit','norm',...
    'position',[0.02 0.29 0.2 0.03],...
    'string','成像面总边长/cm',...
    'fontsize',16);
S.ed_fL = uicontrol('style','edit',...
    'unit','norm',...
    'position',[0.02 0.25 0.1 0.03],...
    'string','1.20000',...
    'fontsize',12);
align([S.tx_fL S.ed_fL],'Center','Fixed',7);


%% 迭代次数
S.tx_XX = uicontrol('style','text',...
    'unit','norm',...
    'position',[0.22 0.035 0.1 0.03],...
    'string','迭代次数 =',...
    'fontsize',16);
S.ed_XX = uicontrol('style','edit',...
    'unit','norm',...
    'position',[0.32 0.035 0.03 0.03],...
    'string','101',...
    'fontsize',12);

%% 显示间隔
S.tx_dXX = uicontrol('style','text',...
    'unit','norm',...
    'position',[0.35 0.035 0.1 0.03],...
    'string','显示间隔 =',...
    'fontsize',16);
S.ed_dXX = uicontrol('style','edit',...
    'unit','norm',...
    'position',[0.45 0.035 0.03 0.03],...
    'string','10',...
    'fontsize',12);

%% 图片文件路径
S.tx_fpath_in = uicontrol('style','text',...
    'unit','norm',...
    'position',[0.2 0.1 0.15 0.03],...
    'string','原图路径：',...
    'fontsize',16);
S.ed_fpath_in = uicontrol('style','edit',...
    'unit','norm',...
    'position',[0.32 0.1 0.45 0.03],...
    'string','.',...
    'fontsize',12);

%% 图片文件名
S.tx_fname_in = uicontrol('style','text',...
    'unit','norm',...
    'position',[0.7 0.1 0.15 0.03],...
    'string','原图名：',...
    'fontsize',16);
S.ed_fname_in = uicontrol('style','edit',...
    'unit','norm',...
    'position',[0.81 0.1 0.1 0.03],...
    'string','heart.png',...
    'fontsize',12);

%% Axis
S.ax3 = axes('units','norm',...
    'position',[0.25 0.2 0.3 0.35]);

S.ax4 = axes('units','norm',...
    'position',[0.6 0.2 0.3 0.35]);

S.ax1 = axes('units','norm',...
    'position',[0.25 0.6 0.3 0.35]);

S.ax2 = axes('units','norm',...
    'position',[0.6 0.6 0.3 0.35]);

%% check 成像面->DOE
S.ch0 = uicontrol('style','check',...
    'units','norm',...
    'position',[0.07 0.4 0.1 0.06],...
    'string','成像面->DOE',...
    'fontsize',12);

%% check 是否负片
S.ch1 = uicontrol('style','check',...
    'units','norm',...
    'position',[0.5 0.02 0.1 0.06],...
    'string','负片',...
    'fontsize',12);
%% check 是否存图
S.ch2 = uicontrol('style','check',...
    'units','norm',...
    'position',[0.7 0.02 0.1 0.06],...
    'string','存图',...
    'fontsize',12);

%% 物面像面换算
S.pb0 = uicontrol('style','push',...
    'units','norm',...
    'position',[0.06 0.35 0.12 0.06],...
    'backgroundcolor',[39/255 177/255 216/255],...
    'HorizontalAlign','Center',...
    'string','物面像面换算',...
    'fontsize',14,'fontweight','bold',...
    'callback',{@pb0_call,S});


%% 载入照片
S.pb1 = uicontrol('style','push',...
    'units','norm',...
    'position',[0.55 0.02 0.12 0.06],...
    'backgroundcolor',[142/255 179/255 238/255],...
    'HorizontalAlign','left',...
    'string','载入照片',...
    'fontsize',14,'fontweight','bold',...
    'callback',{@pb1_call,S});

%% 运算按钮
S.pb2 = uicontrol('style','push',...
    'units','norm',...
    'position',[0.75 0.02 0.15 0.06],...
    'backgroundcolor',[142/255 179/255 238/255],...
    'HorizontalAlign','left',...
    'string','载入并运行',...
    'fontsize',14,'fontweight','bold',...
    'callback',{@pb2_call,S});


function [] = pb0_call(varargin)
% Callback for pushbutton.
[h,S] = varargin{[1,3]};  % Get calling handle and structure.
ch0_in = get(S.ch0,'Value');

wavelength = str2double(get(S.ed_wavelength,'string'));
f = str2double(get(S.ed_f,'string'));
NN = str2double(get(S.ed_NN,'string'));
pixeld = str2double(get(S.ed_pixeld,'string'));
fL = str2double(get(S.ed_fL,'string'));

if(ch0_in==0)
    Changed_fL = (f*1e-2)*(wavelength*1e-9)/(pixeld*1e-6)/(1e-2);%cm
    set(S.ed_fL,'string',num2str(Changed_fL,'%.5f'));
elseif (ch0_in==1)
    Changed_pixeld = (f*1e-2)*(wavelength*1e-9)/(fL*1e-2)/(1e-6);%um
    set(S.ed_pixeld,'string',num2str(Changed_pixeld,'%.5f'));
end







function [] = pb1_call(varargin)
% Callback for pushbutton.
[h,S] = varargin{[1,3]};  % Get calling handle and structure.
fpath_in = get(S.ed_fpath_in,{'string'});  % Get the image path
fpath_in = fpath_in{1};
fname_in  = get(S.ed_fname_in,{'string'}); % Get the image name
fname_in = fname_in{1};
ch1_in = get(S.ch1,'Value');% 负片？
ch2_in = get(S.ch2,'Value');% 存图？
XX = str2double(get(S.ed_XX,'string'));
dXX = str2double(get(S.ed_dXX,'string'));
NN = str2double(get(S.ed_NN,'string'));

%% Pardef

% f面定义
NumLfx = NN;
NumLfy = NN;
% 0面定义
NumL0x = NumLfx;
NumL0y = NumLfy;


%% f面导入图片 NN*NN
fsavepath = [fpath_in,'\',fname_in,'_save'];
fname = fname_in;

heart = imread([fpath_in,'\',fname_in]);
hahaha = size(heart);
heartX = hahaha(1);
heartY = hahaha(2);
if heartX>=heartY
    heart= imresize(heart,[NN floor(heartY/heartX*NN)]);
else
    heart= imresize(heart,[floor(heartX/heartY*NN) NN]);
end

heart = sum(heart,3);
heart = heart/max(max(heart));
[heartX,heartY] = size(heart);
if ch1_in == 1
    heart = 1-heart;
end

TargetAmplitude = zeros(NumLfx,NumLfy);
if(heartX>=heartY)
    TargetAmplitude(1:heartX,(1:heartY)+floor((heartX-heartY)/2)) = heart;
else
    TargetAmplitude((1:heartX)+floor((heartY-heartX)/2),1:heartY) = heart;
end



%% 0面定义
Pattern0 = zeros(NumL0x,NumL0y);
Source0 = ones(NumL0x,NumL0y);


axes(S.ax1);
cla
hold on
title(['Source0']);
imagesc(Source0);
cmap = gray(256);
colormap(gray);
colorbar();
axis image

axes(S.ax2);
cla

axes(S.ax3);
cla


axes(S.ax4);
cla;
hold on
title(['TargetAmplitude']);
imagesc(flip(rot90(rot90(TargetAmplitude)),2));
cmap = gray(256);
colormap(gray);
colorbar();
axis image
xlim([0,NN])
ylim([0,NN])










function [] = pb2_call(varargin)
% Callback for pushbutton.
[h,S] = varargin{[1,3]};  % Get calling handle and structure.
fpath_in = get(S.ed_fpath_in,{'string'});  % Get the image path
fpath_in = fpath_in{1};
fname_in  = get(S.ed_fname_in,{'string'}); % Get the image name
fname_in = fname_in{1};
ch1_in = get(S.ch1,'Value');% 负片？
ch2_in = get(S.ch2,'Value');% 存图？
XX = str2double(get(S.ed_XX,'string'));
dXX = str2double(get(S.ed_dXX,'string'));
NN = str2double(get(S.ed_NN,'string'));

%% Pardef
% f面定义
NumLfx = NN;
NumLfy = NN;
% 0面定义
NumL0x = NumLfx;
NumL0y = NumLfy;


%% f面导入图片 NN*NN
fsavepath = [fpath_in,'\',fname_in,'_save'];
if(ch2_in==1)
    mkdir(fsavepath);
end
fname = fname_in;
%     fpath = [pwd,'\1_HeartImage'];
%     fname = ['heart.png'];

heart = imread([fpath_in,'\',fname_in]);
hahaha = size(heart);
heartX = hahaha(1);
heartY = hahaha(2);
if heartX>=heartY
    heart= imresize(heart,[NN floor(heartY/heartX*NN)]);
else
    heart= imresize(heart,[floor(heartX/heartY*NN) NN]);
end

heart = sum(heart,3);
heart = heart/max(max(heart));
[heartX,heartY] = size(heart);
if ch1_in == 1
    heart = 1-heart;
end



TargetAmplitude = zeros(NumLfx,NumLfy);

if(heartX>=heartY)
    TargetAmplitude(1:heartX,(1:heartY)+floor((heartX-heartY)/2)) = heart;
else
    TargetAmplitude((1:heartX)+floor((heartY-heartX)/2),1:heartY) = heart;
end


%% 0面定义
Pattern0 = zeros(NumL0x,NumL0y);
Source0 = ones(NumL0x,NumL0y);

%% GS Agorithm
TargetPhase = pi*(rand(NumLfx,NumLfy)-0.5)*2; %初始相位
DOEphase = angle(ifft2(fftshift(TargetAmplitude).*exp(1i*TargetPhase)));
DOEphase = floor((DOEphase+pi)/(2*pi)*256)/256;
for kk =1:XX
    kk;
    AA = fft2(Source0.*exp(1i*DOEphase));
    TargetPhase = angle(AA);
    
    if (mod(kk-1,dXX)==0)
        kk
        %                     figure(888);
        
        %                     subplot(2,2,1);
        axes(S.ax1);
        cla
        hold on
        title(['Source0']);
        imagesc(Source0);
        cmap = gray(256);
        colormap(gray);
        colorbar();
        axis image
        
        %                     subplot(2,2,2);
        axes(S.ax2);
        cla
        hold on
        title(['IterationTimes=',num2str(kk),'---DOEphase']);
        imagesc(DOEphase);
        cmap = gray(256);
        colormap(gray);
        colorbar();
        axis image
        
        %                     subplot(2,2,3);
        axes(S.ax3);
        cla
        hold on
        title(['IterationTimes=',num2str(kk),'---abs(AA)']);
        imagesc(flip(rot90(rot90(abs(fftshift(AA)))),2));
        cmap = gray(256);
        colormap(gray);
        colorbar();
        axis image
        xlim([0,NN])
        ylim([0,NN])
        
        %                     subplot(2,2,4);
        axes(S.ax4);
        cla;
        hold on
        title(['TargetAmplitude']);
        imagesc(flip(rot90(rot90(TargetAmplitude)),2));
        cmap = gray(256);
        colormap(gray);
        colorbar();
        axis image
        xlim([0,NN])
        ylim([0,NN])
        
        
        if(ch2_in==1)
            print(gcf,'-dpng',[fsavepath,'\',fname,'IterationTimes',num2str(kk),'.png'])
        end
        
        %         print(gcf,'-dpng',[fsavepath,'\',fname,'IterationTimes',num2str(kk),'.png'])
        %         close(gcf);
        
    end
    
    DOEphase = angle(ifft2(fftshift(TargetAmplitude).*exp(1i*TargetPhase)));
    DOEphase = floor((DOEphase+pi)/(2*pi)*256)/256*2*pi;
end


if(ch2_in==1)
    save([fsavepath,'\_',fname,'IterationTimes',num2str(kk),'---',num2str(NumLfx),...
        '-',num2str(NumLfy),'DOEPhase.mat'],'DOEphase')
    %% 直接写入图
    imwrite([flip(rot90(rot90(flip(TargetAmplitude,1)/max(max(TargetAmplitude)))),2),...
        DOEphase/max(max(DOEphase))],...
        [fsavepath,'\_',fname,'IterationTimes',num2str(kk),'---',num2str(NumLfx),...
        '-',num2str(NumLfy),'DOEPhase.png'],'PNG');
end