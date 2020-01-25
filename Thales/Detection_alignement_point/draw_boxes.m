function draw_boxes(align_m)

Nm = size(align_m,2);
% draw
if (Nm ~= 0)
    % alignments
    for i=0:Nm-1
        x1  = align_m(8*i + 1);
        y1  = align_m(8*i + 2);
        x2  = align_m(8*i + 3);
        y2  = align_m(8*i + 4);
        w   = align_m(8*i + 5);
        l   = align_m(8*i + 6);
        b   = align_m(8*i + 7);
        nfa = align_m(8*i + 8);
        dx = x2-x1;
        dy = y2-y1;
        len = sqrt( dx*dx + dy*dy );
        dx = dx/len;
        dy = dy/len;
        
        hw = w / 2.0; % half width
        hl = l / 2.0; % half local window width
        
        % alignment rectangle
        line([x1-dy*hw x2-dy*hw],[y1+dx*hw y2+dx*hw],'Color','green','LineWidth',1);
        line([x1+dy*hw x2+dy*hw],[y1-dx*hw y2-dx*hw],'Color','green','LineWidth',1);
        line([x1-dy*hw x1+dy*hw],[y1+dx*hw y1-dx*hw],'Color','green','LineWidth',1);
        line([x2-dy*hw x2+dy*hw],[y2+dx*hw y2-dx*hw],'Color','green','LineWidth',1);
        
        % local window
        if( l > 0.0 )
            line([x1-dy*hl x2-dy*hl],[y1+dx*hl y2+dx*hl],'Color','blue','LineWidth',2);
            line([x1+dy*hl x2+dy*hl],[y1-dx*hl y2-dx*hl],'Color','blue','LineWidth',2);
            line([x1-dy*hl x1+dy*hl],[y1+dx*hl y1-dx*hl],'Color','blue','LineWidth',2);
            line([x2-dy*hl x2+dy*hl],[y2+dx*hl y2-dx*hl],'Color','blue','LineWidth',2);
        end
        
        % boxes
        if( b > 0.0 )
            hb = len / (b+1) / 2.0; % half box side
            for j=1:b
                line([x1-dy*hw+dx*hb*(1+2*j) x1+dy*hw+dx*hb*(1+2*j)],[y1+dx*hw+dy*hb*(1+2*j) y1-dx*hw+dy*hb*(1+2*j)], 'Color','blue','LineWidth',1);
            end
        end   
    end

end