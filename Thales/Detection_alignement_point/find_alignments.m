
% input:  point     - input points (including extended points at the end)
% N         - number of points
% Next      - number of extended points

% output: pointer to *Nout 8-tuples
% x1 y1 x2 y2 width local_window_width num_boxes log_nfa

function align = find_alignments(point, N, Next, X, Y)
align = [];
n_align = 0;
max_w_over_len = 10;
n_w = 8;
n_l = 8;
% number of tests
logNT = log10(N) + log10(N-1) - log10(2) + log10(n_w) + log10(n_l) + log10(N)/2.0;

% loop over all pair of points
for i = 1:N-1
    for j = i+1:N
        x1 = point(i,1);	y1 = point(i,2);
        x2 = point(j,1);	y2 = point(j,2);
        dx = x2-x1;		dy = y2-y1;
        len = sqrt(dx*dx+dy*dy);
        if( len == 0 )
            continue;
        end
        dx = dx/len ;		dy = dy/len;
        
        % loop over width of our alignement
        w = len/max_w_over_len;
        for ww = 1:n_w
            % loop over size of our rectangle = local window
            l = len;
            for ll = 1:n_l
                
                % check if local window goes out of domain.
                % to speed-up, extended points are only used when the local window exceed the domain.
                if( (x1 - dy*l/2) < 0 || (x1 - dy*l/2) >= X || (x1 + dy*l/2) < 0 || (x1 + dy*l/2) >= X ||...
                        (x2 - dy*l/2) < 0 || (x2 - dy*l/2) >= X || (x2 + dy*l/2) < 0 || (x2 + dy*l/2) >= X ||...
                        (y1 + dx*l/2) < 0 || (y1 + dx*l/2) >= Y || (y1 - dx*l/2) < 0 || (y1 - dx*l/2) >= Y ||...
                        (y2 + dx*l/2) < 0 || (y2 + dx*l/2) >= Y || (y2 - dx*l/2) < 0 || (y2 - dx*l/2) >= Y )
                    NN = Next;
                else
                    NN = N;
                end
                
                % count number of points in the strips
                n1 = 0; n2 = 0; nn = 0;
                for m = 1:NN
                    if( (m==i) || (m==j) )
                        continue;
                    end % do not count i and j : the 2 points that defines the alignment
                    
                    % compute coordinates relative to alignments
                    xx =  dx * (point(m,1)-x1) + dy * (point(m,2)-y1);
                    yy = -dy * (point(m,1)-x1) + dx * (point(m,2)-y1);
                    
                    % count local window points
                    if ( (xx < 0) || (xx >= len) ) % out of the rectangle in x
                        continue; % pass to next iteration ...
                    end
                    if ( (yy < -l/2) || (yy > l/2) ) % out of the rectangle in y
                        continue;
                    end
                    if ( yy < -w/2 ) % the bottom of the rectangle
                        n1 = n1+1;
                    elseif ( yy > w/2 ) % the top of the rectangle
                        n2 = n2+1;
                    elseif ( m < N ) % middle of the rectangle and not out of domain
                        nn = nn+1;
                        in_align(nn) = xx;
                    end
                end
                
                % number of boxes in the alignment: from 1 to N-2
                for k = max(floor(nn/2),1) : min(2*nn,N-2)
                    % compute box side
                    box = len/(k+1);
                    % free used boxes
                    used(1:k) = 0;
                    % count point in boxes
                    for m = 1:nn
                        b = floor((in_align(m) - box/2)/box);
                        if ( b >= 0 && b < k )
                            used(b+1) = used(b+1)+1;
                        end
                    end
                    
                    % compute NFA
                    n = 2 * max(n1,n2) + nn ;
                    p0 = w * box / len / l ;
                    p1 = 1 - (1 - p0)^ n ;
                    n_box = 0 ;
                    for v = 1:k
                        if( used(v) ~= 0 )
                            n_box = n_box+1;
                        end
                    end
                    logNFA = logNT + log_bin(k, n_box, p1);
                 
                    % store result
                    if (logNFA < 0)
                         align(8*n_align+1) = x1;
                         align(8*n_align+2) = y1;
                         align(8*n_align+3) = x2;
                         align(8*n_align+4) = y2;
                         align(8*n_align+5) = w;
                         align(8*n_align+6) = l;
                         align(8*n_align+7) = k;
                         align(8*n_align+8) = logNFA;
                         n_align = n_align + 1 ;
                    end
                    
                end
                l = l/sqrt(2);
            end
            w = w/sqrt(2);
        end 
        
    end
    
end

end









