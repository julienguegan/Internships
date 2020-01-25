clear all
clf

cpt = 1;
% data rate :
% FT (false transient) : 80000/s donc 80/ms
% TT en cluster de 2x2 : 4000/s donc 4/ms

% x = 1:10:10;
% for N = 1:10:100
N = 0; % nombre de points aleatoires

% s = rng; % pour garder la meme distribution aleatoire quand on a besoin de boucler sur les points pour les tests
% point(1:N,1) = round(rand(1,N)*1000+1); % coord en x
% point(1:N,2) = round(rand(1,N)*1170+1); % coord en y
% rng(s);

% circle
for k = 1:2
    xC = (k-1)*500+200; yC = (k-1)*500+200; R = k*40; s = round(R*1.5);
    th = rand(1,s)*2*pi; r = rand(1,s)*R;
    x = round(r.*cos(th)+xC); y = round(r.*sin(th)+yC); 
    xpos = x((x>0)+(y>0)+(x<1000)+(y<1170)==4); ypos = y((x>0)+(y>0)+(x<1000)+(y<1170)==4); 
    s = length(xpos);
    point(N+1:N+s,1) = xpos; % coord en x
    point(N+1:N+s,2) = ypos; % coord en y
    N = N+s;
end

box = 'box';
algo = 'A Contrario';

% setup alignment
lp1 = 0;
debut = N+1;
lp0 = 15;
point(debut:debut-1+lp0,1) = linspace(200,400,15);
point(debut:debut-1+lp0,2) = linspace(400,800,15);
csmprtPoint(1:lp0,1) = linspace(200,300,15);
csmprtPoint(1:lp0,2) = linspace(300,800,15);
point(debut+lp0:debut+2*lp0-1,1) = linspace(200,700,15);
point(debut+lp0:debut+2*lp0-1,2) = linspace(600,300,15);
csmprtPoint(lp0+1:2*lp0,1) = linspace(200,300,15);
csmprtPoint(lp0+1:2*lp0,2) = linspace(600,300,15);

tic;
[pointInAlign_AContrario,align_m] = point_alignments(point);
timeAContrario(cpt) = toc;
% -------- check output --------- %
if ( size(pointInAlign_AContrario,1) > 0 )
    accuracyAContrario(cpt) = size(intersect(csmprtPoint,pointInAlign_AContrario,'rows'),1)/size(csmprtPoint,1)*100;
end


% -------- draw --------- %
figure(cpt);
hold on
scatter(point(:,1),point(:,2),30,'k','filled');
if (size(pointInAlign_AContrario,1)>0)
    if (strcmp(algo,'A Contrario') && strcmp(box,'box'))
        draw_boxes(align_m);
    end
    scatter(pointInAlign_AContrario(:,1),pointInAlign_AContrario(:,2),30,'r','filled');
    title([algo ' - ' num2str(timeAContrario(cpt)) ' s - ' num2str(accuracyAContrario(cpt)) ' %']);
end
hold off


cpt = cpt+1;



