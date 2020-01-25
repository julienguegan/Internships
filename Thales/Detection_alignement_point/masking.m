% input:  point     - input points (including extended points at the end)
% N         - number of points
% Next      - number of extended points
% align     - list of alignments
% Na        - number of input alignments

% output: pointer to *Nout 8-tuples
% x1 y1 x2 y2 width local_window_width num_boxes log_nfa

function align_m = masking(point, N, Next, align, Na, X, Y)
align_m = [];
n_w = 8;
n_l = 8;
% number of tests
logNT = log10(N) + log10(N-1) - log10(2) + log10(n_w) + log10(n_l) + log10(N)/2;

% sort alignments by NFA
align = reshape(align,[8,Na]);
align = sortrows(align', 8);
align = reshape(align',[],1);

n_out = 0;
% first alignment is directly added to the output
align_m(1,1) = align(1);
align_m(2,1) = align(2);
align_m(3,1) = align(3);
align_m(4,1) = align(4);
align_m(5,1) = align(5);
align_m(6,1) = align(6);
align_m(7,1) = align(7);
align_m(8,1) = align(8);
n_out = n_out+1;

% loop over ordered alignments
for i = 1:(Na-1)
    
    % alignment data
    x1 = align(8*i+1);
    y1 = align(8*i+2);
    x2 = align(8*i+3);
    y2 = align(8*i+4);
    w  = align(8*i+5);
    l  = align(8*i+6);
    k  = align(8*i+7);
    dx = x2-x1;
    dy = y2-y1;
    len = sqrt(dx*dx+dy*dy);
    dx = dx/len;
    dy = dy/len;
    box = len / (k+1);
    
    % try masking the alignment by already detected ones
    max_logNFA = align(8*i+8);
    for j = 1:n_out
        % masking candidate data
        m_x1 = align_m(1,j);
        m_y1 = align_m(2,j);
        m_x2 = align_m(3,j);
        m_y2 = align_m(4,j);
        m_w  = align_m(5,j);
        m_l  = align_m(6,j);
        m_dx = m_x2-m_x1;
        m_dy = m_y2-m_y1;
        m_len = sqrt(m_dx*m_dx+m_dy*m_dy);
        m_dx = m_dx/m_len;
        m_dy = m_dy/m_len;
        % free used boxes
        used(1:k) = 0;
        
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
        
        % count points
        n1 = 0; n2 = 0; n3 = 0;
        for m = 1:NN
            % do not count defining points
            if( (point(m,1) == x1) && (point(m,2) == y1) )
                continue;
            end
            if( (point(m,1) == x2) && (point(m,2) == y2) )
                continue;
            end
            
            % compute coordinates relative to alingments
            xx =  dx * (point(m,1)-x1) + dy * (point(m,2)-y1);
            yy = -dy * (point(m,1)-x1) + dx * (point(m,2)-y1);
            
            % count points in the local window
            if( (xx < 0) || (xx >= len) )
                continue;
            end
            if( (yy < -l/2) || (yy > l/2) )
                continue;
            end
            if( yy < -w/2 )
                n1 = n1+1;
                continue; % point in local window but not in alignment
            end
            if( yy >  w/2 )
                n2 = n2+1;
                continue; % point in local window but not in alignment
            end
            if( m >= N )
                continue;
            end % do not count ext points
            n3 = n3+1;
            
            % masked if the point belong to masking candidate
            m_xx =  m_dx * (point(m,1)-m_x1) + m_dy * (point(m,2)-m_y1);
            m_yy = -m_dy * (point(m,1)-m_x1) + m_dx * (point(m,2)-m_y1);
            if( (abs(m_yy) < m_w/2) && (m_xx>=0) && (m_xx<m_len) )
                continue;
            end
            
            % count if the points that belong to one box
            b = floor( (xx - box/2) / box );
            if( (b>=0) && (b<k) )
                used(b+1) = used(b+1)+1;
            end
        end
        
        % compute NFA after masking
        n = 2 * max(n1,n2) + n3;
        p0 = w * box / len / l;
        p1 = 1 - (1 - p0)^n ;
        n_box = 0;
        for v = 1:k
            if( used(v) ~= 0 )
                n_box = n_box+1;
            end
        end
        logNFA = logNT + log_bin( k, n_box, p1 );
        if( logNFA > max_logNFA )
            max_logNFA = logNFA;
        end
    end
    
    if( max_logNFA > 0 )
        continue;
    end
    
    % store alignment to output
    n_out = n_out+1;
    align_m(1,n_out) = align(8*i+1);
    align_m(2,n_out) = align(8*i+2);
    align_m(3,n_out) = align(8*i+3);
    align_m(4,n_out) = align(8*i+4);
    align_m(5,n_out) = align(8*i+5);
    align_m(6,n_out) = align(8*i+6);
    align_m(7,n_out) = align(8*i+7);
    align_m(8,n_out) = align(8*i+8);
    
    
end

end
