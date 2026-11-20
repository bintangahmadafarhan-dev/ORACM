function Web    
    clear; close all; clc;
    list = dir('images\'); list = list(4:end,1);
    
    hFig = figure;
    filename = 'animation.gif';
    
    img = imread(['images\',list(1,1).name]);
    if size(img,3)>1,img = rgb2gray(img);end
    img = double(img);
    
    Img = img;
    [row,col] = size(Img);
    u = -ones(row,col);
    u(10:row-10,10:col-10) = 1;   
    sigma = 1; G = fspecial('gaussian', 5, sigma);
    delt = 1; mu = 10;itr = 0; 
    objPos = u >= 0; objNeg = ~objPos;
    Area1 = sum(u(:)); Area2 =0;
    
    subplot(3, 1, 1);
    imshow(Img,[]); hold on; 
    contour(u, [0 0], 'r','LineWidth',4);hold on;
    contour(u, [0 0], 'g','LineWidth',1.3);hold on;
    title(['ACMwithSBGFRLS : ',int2str(itr),'th iteration']);axis off; 
        
    frame = getframe(hFig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);

    tic    
    while abs(Area2-Area1)>0      
        c1 = sum(sum(Img.*objNeg))/sum(sum(objNeg));
        c2 = sum(sum(Img.*objPos))/sum(sum(objPos));     
        spf = Img - (c1 + c2)/2;
        spf = spf/(max(abs(spf(:))));    
        [ux, uy] = gradient(u);
        u = u + delt*(mu*spf.*sqrt(ux.^2 + uy.^2));
        u = (u >= 0) - ( u< 0);
        u = conv2(u, G, 'same'); 
        objPos = u >= 0; objNeg = ~objPos;  
        Area1= Area2;  Area2 =sum(objPos(:)); 
        itr = itr + 1;
        subplot(3, 1, 1);
        imshow(Img,[]); hold on; 
        contour(u, [0 0], 'r','LineWidth',4);hold on;
        contour(u, [0 0], 'g','LineWidth',1.3);hold on;
        title(['ACMwithSBGFRLS : ',int2str(itr),'th iteration']);axis off; 
        drawnow; 
    
        frame = getframe(hFig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);        
    end  

    subplot(3, 1, 2);
    I = img;
    k = 10; [row,col] = size(I);
    mask = zeros(row,col);
    mask (k:row-k,k:col-k)=1;    
    mu=0.2; itr = 0;
    phi0 = bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5;         
    objPos = phi0 >= 0; objNeg = ~objPos;
    Area1 = sum(objNeg(:)); Area2 =sum(objPos(:));

    imshow(I);         
    contour(flipud(phi0), [0 0], 'r','LineWidth',4); hold on; 
    contour(flipud(phi0), [0 0], 'g','LineWidth',1.3); hold on; 
    title(['C-V: ',int2str(itr),'th iteration']);axis off;         
    drawnow;

    frame = getframe(hFig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);      

    tic
    while abs(Area2-Area1)>0.3   
        inidx = find(phi0>=0); % frontground index
        outidx = find(phi0<0); % background index
        force_image = 0; % initial image force for each layer                   
        c1 = sum(sum(I.*Heaviside(phi0)))/(length(inidx)+eps); % average inside of Phi0
        c2 = sum(sum(I.*(1-Heaviside(phi0))))/(length(outidx)+eps); % verage outside of Phi0
        force_image=force_image-(I-c1).^2+(I-c2).^2; 
        kphi0 = kappa(phi0); maxkphi0 = max(max(abs(kphi0)));
        force = force_image + mu*kphi0./maxkphi0;
        force = force./(max(max(abs(force))));        
        phi0 = phi0+force;        
        objPos=phi0(abs(phi0)<=.5);       
        Area1= Area2;  Area2 =sum(objPos(:));        
        itr = itr + 1; 
        subplot(3, 1, 2);          
        imshow(I,'initialmagnification','fit','displayrange',[0 255]);hold on;
        contour(phi0, [0 0], 'r','LineWidth',4);hold on; 
        contour(phi0, [0 0], 'g','LineWidth',1.3);hold on; 
        title(['C-V : ',int2str(itr),'th iteration']);axis off;         
        drawnow;

        frame = getframe(hFig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);            
    end

    subplot(3, 1, 3);
    Img = img;
    warning('off'); %#ok<WNOFF>
    [row,col] = size(Img);
    phi = ones(row,col); k = 10; 
    phi(k:row-k,k:col-k) = -1;
    obj = - phi; objPos = obj >= 0; objNeg = ~objPos;
    itr = 0; 
    Area1 = sum(objNeg(:)); Area2 =sum(objPos(:));

    imshow(Img,[]); hold on; 
    contour(obj, [0 0], 'r','LineWidth',4);hold on; 
    contour(obj, [0 0], 'g','LineWidth',1.3);hold on; 
    title(['Proposed Method (ORACM):',int2str(itr),'th iteration']);axis off; 

    frame = getframe(hFig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);       

    tic    
    while abs(Area2-Area1)>0                                            
        c1 = sum(sum(Img.*objNeg))/sum(sum(objNeg));
        c2 = sum(sum(Img.*objPos))/sum(sum(objPos)); 
        nImg = Img - (c1 + c2)/2;
        obj = nImg /max(abs(nImg(:)));            
        objPos = obj >= 0; objNeg = ~objPos;
        obj = objPos - objNeg;            
        Area1= Area2;  Area2 =sum(objPos(:));    
        itr = itr + 1;

        subplot(3, 1, 3);
        imshow(Img,[]); hold on; 
        contour(obj, [0 0], 'r','LineWidth',4);hold on; 
        contour(obj, [0 0], 'g','LineWidth',1.3);hold on; 
        title(['Proposed Method (ORACM):',int2str(itr),'th iteration']);axis off; 
        drawnow;

        frame = getframe(hFig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);           
    end 
    imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 2);
    close(hFig);
end

function H=Heaviside(z)
    Epsilon=10^(-5);
    H=zeros(size(z,1),size(z,2));
    idx1=find(z>=Epsilon);
    idx2=find(z<Epsilon & z>-Epsilon);
    H(idx1)=1;
    for i=1:length(idx2)
        H(idx2(i))=1/2*(1+z(idx2(i))/Epsilon+1/pi*sin(pi*z(idx2(i))/Epsilon));
    end;
end

function KG = kappa(I)  
    I = double(I);
    [m,n] = size(I);
    P = padarray(I,[1,1],1,'pre');
    P = padarray(P,[1,1],1,'post');

    % central difference
    fy = P(3:end,2:n+1)-P(1:m,2:n+1);
    fx = P(2:m+1,3:end)-P(2:m+1,1:n);
    fyy = P(3:end,2:n+1)+P(1:m,2:n+1)-2*I;
    fxx = P(2:m+1,3:end)+P(2:m+1,1:n)-2*I;
    fxy = 0.25.*(P(3:end,3:end)-P(1:m,3:end)+P(3:end,1:n)-P(1:m,1:n));
    G = (fx.^2+fy.^2).^(0.5);
    K = (fxx.*fy.^2-2*fxy.*fx.*fy+fyy.*fx.^2)./((fx.^2+fy.^2+eps).^(1.5));
    KG = K.*G;
    KG(1,:) = eps;
    KG(end,:) = eps;
    KG(:,1) = eps;
    KG(:,end) = eps;
    KG = KG./max(max(abs(KG)));
end