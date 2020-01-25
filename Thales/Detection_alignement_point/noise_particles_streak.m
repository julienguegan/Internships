function point = noise_particles_streak(row, column,length_min, thick,  angles_num)
% Cosmic Particle noise


thick_max=max(thick);

angles_vect=rand([1,angles_num])*180;

length_vect=length_min+round(rand([1,angles_num])/10*length_min);

x_array=zeros(4,1);
y_array=zeros(4,1);
point_x = [];
point_y = [];
for cycle1=1:angles_num
    
    rot_angle=angles_vect(cycle1)*pi/180;  % deg2rad in R2015b deg2rad(angles_vect(cycle1))
    
    height=round(length_vect(cycle1)*sin(rot_angle));
    
    ok=false;
    
    while ok==false
        
        x1=row-round(rand(1)*(row-1));
        y1=round(rand(1)*(column-1))+1;
        
        x2=x1-height;
        y2=round(y1+cos(rot_angle)*length_vect(cycle1));
        
        
        if (x2>=1+thick_max) && (x2<=row-thick_max) && (y2>=1+thick_max) && (y2<=column-thick_max)
            ok=true;
        end
    end
    
    y_array(1)=y1+sin(rot_angle)*thick(1)/2;
    y_array(2)=y2+sin(rot_angle)*thick(2)/2;
    y_array(3)=y2-sin(rot_angle)*thick(2)/2;
    y_array(4)=y1-sin(rot_angle)*thick(1)/2;
    % x_array et y_array = bord du polygone
    x_array(1)=x1+cos(rot_angle)*thick(1)/2;
    x_array(2)=x2+cos(rot_angle)*thick(2)/2;
    x_array(3)=x2-cos(rot_angle)*thick(2)/2;
    x_array(4)=x1-cos(rot_angle)*thick(1)/2;
    
    [xq,yq] = meshgrid((1:row)',(1:column)');
    % on teste tous les couples dans [1:1000,1:1170]
    in = inpolygon(xq,yq,x_array,y_array); % dit quel point du plan sont dans le polygone
    
    point_x = [point_x xq(in)]; % concatene if we ask several alignment
    point_y = [point_y yq(in)];
end

point(:,1) = point_x;
point(:,2) = point_y;

end

 % rand(500,1)*5 : 500 doubles entre 0 et 5



