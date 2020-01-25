clear all
clf

cpt = 1;
% data rate :
% FT (false transient) : 80000/s donc 80/ms
% TT en cluster de 2x2 : 4000/s donc 4/ms

% x = 1:10:10;
% for N = 1:10:100
N = 30; % nombre de points aleatoires

% s = rng; % pour garder la meme distribution aleatoire quand on a besoin de boucler sur les points pour les tests
point(1:N,1) = round(rand(1,N)*1000+1); % coord en x
point(1:N,2) = round(rand(1,N)*1170+1); % coord en y
% rng(s);
for k = 1:3
    xC = rand()*1000+1; yC = rand()*1170+1; R = rand()*40+20; s = round(R/1.5);
    th = rand(1,s)*2*pi; r = rand(1,s)*R;
    x = round(r.*cos(th)+xC); y = round(r.*sin(th)+yC); 
    xpos = x((x>0)+(y>0)+(x<1000)+(y<1170)==4); ypos = y((x>0)+(y>0)+(x<1000)+(y<1170)==4); 
    s = length(xpos);
    point(N+1:N+s,1) = xpos; % coord en x
    point(N+1:N+s,2) = ypos; % coord en y
    N = N+s;
end


type = 'particule';
box = 'nobox';
algo = 'A Contrario';

% setup particle noise
cosmic_prt.angles_num = 1;
cosmic_prt.thick = [2 2]; % [@start @end]
cosmic_prt.length = round(rand()*20+30);
cosmic_prt.number = round(rand()*5+5); % entre 5 et 10 particule
lp1 = 0;
debut = N+1;
for i = 1:cosmic_prt.number
    particules = noise_particles_streak(1000, 1170, cosmic_prt.length, cosmic_prt.thick, cosmic_prt.angles_num) ;
    lp0 = length(particules(:,1));
    point(debut:debut-1+lp0,1) = particules(:,1);
    point(debut:debut-1+lp0,2) = particules(:,2);
    csmprtPoint(1+lp1*(i-1):lp0+lp1*(i-1),1) = particules(:,1);
    csmprtPoint(1+lp1*(i-1):lp0+lp1*(i-1),2) = particules(:,2);
    lp1 = lp0;
    debut = debut+lp1;
end

tic;
[pointInAlign_AContrario,align_m] = point_alignments_for_csmprt(point,cosmic_prt.length-1,max(cosmic_prt.thick));
timeAContrario(cpt) = toc;
% -------- check output --------- %
if ( size(pointInAlign_AContrario,1) > 0 )
    accuracyAContrario(cpt) = size(intersect(csmprtPoint,pointInAlign_AContrario,'rows'),1)/size(csmprtPoint,1)*100;
end


% -------- draw --------- %
figure(cpt);
%clf(3);
subplot(1,2,1)
hold on
scatter(point(:,1),point(:,2),12,'k','filled');
if (size(pointInAlign_AContrario,1)>0)
    if (strcmp(algo,'A Contrario') && strcmp(box,'box'))
        draw_boxes(align_m);
    end
    scatter(pointInAlign_AContrario(:,1),pointInAlign_AContrario(:,2),12,'r','filled');
end

title([algo ' - ' num2str(timeAContrario(cpt)) ' s - ' num2str(accuracyAContrario(cpt)) ' %']);
hold off


tic;
pointInAlign_Hough = leonardo_test(point,cosmic_prt.length-1,max(cosmic_prt.thick)+1,cosmic_prt.number+1);
timeHough(cpt) = toc;
% -------- check output --------- %
if ( size(pointInAlign_Hough,1) > 0 )
    accuracyHough(cpt) = size(intersect(csmprtPoint,pointInAlign_Hough,'rows'),1)/size(csmprtPoint,1)*100;
end
% -------- draw --------- %
subplot(1,2,2)
hold on
scatter(point(:,1),point(:,2),12,'k','filled');
if (size(pointInAlign_Hough,1)>0)
    scatter(pointInAlign_Hough(:,1),pointInAlign_Hough(:,2),12,'r','filled');
end
title(['Hough - ' num2str(timeHough(cpt)) ' s - ' num2str(accuracyHough(cpt)) ' %']);
hold off

cpt = cpt+1;

% end
%

% h=figure();
% clf();
% ax1 = subplot(1,2,1);
% hold on
% plot(x,accuracyHough,'.r',x, accuracyAContrario,'.b');
% pH = polyfit(x,accuracyHough,1); pAC = polyfit(x,accuracyAContrario,1);
% yH = polyval(pH,x); yAC = polyval(pAC,x);
% plot(x,yH,'-r'); plot(x,yAC,'-b');
% hold off
% title('Accuracy');
% xlabel('nombre de points bruits');
% legend('Hough','A Contrario');
%
%ax2 = subplot(1,2,2);
% hold on
% plot(x,timeHough,'.r', x, timeAContrario,'.b');
% pH = polyfit(x,timeHough,1); pAC = polyfit(x,timeAContrario,2);
% yH = polyval(pH,x); yAC = polyval(pAC,x);
% plot(x,yH,'-r'); plot(x,yAC,'-b');
% hold off
% title('Time');
% xlabel('nombre de points bruits');
% legend('Hough','A Contrario');


