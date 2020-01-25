% input:
%       point     - input points (including extended points at the end)
%       N         - number of points
%       X,Y       - domain size [0,X] x (0,Y]

% output:
%       pointer to *Nout 8-tuples
%       x1 y1 x2 y2 width local_window_width num_boxes log_nfa
function [point_in_align,align_m] = point_alignments_for_csmprt(point,length_min,thick_max)
point_in_align = [];
align_m = [];
% domain
X = 1000;
Y = 1170;

N = size(point,1);

Xsign =   [ 1, -1,  1, -1, -1, -1, -1,  1, -1 ];
Ysign =   [ 1, -1, -1, -1,  1,  1, -1, -1, -1 ];
Xoffset = [ 0,  0,  0,  2,  0,  2,  0,  0,  2 ];
Yoffset = [ 0,  0,  0,  0,  0,  0,  2,  2,  2 ];

% extended points to handle domain border: symmtric extension
Next = 9*N;
for i = 1:9
    for j = 1:N
        point_ext((i-1)*N + j,1) = Xsign(i) * point(j,1) + X*Xoffset(i);
        point_ext((i-1)*N + j,2) = Ysign(i) * point(j,2) + Y*Yoffset(i);
    end
end

% find all meaningful alignments
align = find_alignments_for_csmprt(point_ext,N,Next,X,Y,length_min,thick_max);
Na = length(align)/8;
disp(['find_alignments : ', num2str(Na), '  t : ', num2str(toc),' s']);

align = reshape(align,[8,Na]);
align = sortrows(align', [-6 8]);
[C,ia,ic] = unique(align(:,1:4),'rows','first');
align = align(ia,:);
align = reshape(align',[],1);
Na = length(align)/8;

% reduce redundancy
align_m = [];
if( Na>=1 )
    align_m = masking_for_csmprt(point_ext,N,Next,align,Na,X,Y);
end
    
Nm = size(align_m,2);
disp(['masking : ', num2str(Nm), '  t : ', num2str(toc),' s']);
if (Nm>0)
    nn = 0;
    for i = 1:Nm % loop over the different alignment
        x1 = align_m(1,i); x2 = align_m(3,i); y1 = align_m(2,i); y2 = align_m(4,i);
        dx = x2-x1; dy = y2-y1;
        len = sqrt(dx*dx+dy*dy);
        dx = dx/len; dy = dy/len;
        w = align_m(5,i);
        for m = 1:N  % find the point that are part of the alignment
            xx =  dx * (point(m,1)-x1) + dy * (point(m,2)-y1);
            yy = -dy * (point(m,1)-x1) + dx * (point(m,2)-y1);
            if (xx >= 0 && xx <= len && -w/2 <= yy && yy <= w/2)  % middle of the rectangle
                nn = nn+1;
                point_in_align(nn,1) = point(m,1); % can store several time the same point
                point_in_align(nn,2) = point(m,2);
            end
        end
        
    end
    point_in_align = unique(point_in_align,'rows');
end


end