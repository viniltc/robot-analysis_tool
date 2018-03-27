       %get and save the current size and units
        origFunits = get(gcf,'units');
        set(gcf,'units','inch');
        fpos = get(gcf,'position');
        origPpos = get(gcf,'paperposition');
        origPunits = get(gcf,'paperunits');
        set(gcf,'paperunits','inch');
        psize = get(gcf, 'papersize');

        %does the graph fit on paper
        if fpos(3) > psize(1) | fpos(4) > psize(2) 
            disp( 'WARNING: Graph to large to print in WYSIWYG mode.' )
        else
            llx = (psize(1) - fpos(3)) /2;
            lly = (psize(2) - fpos(4)) /2;
            set( gcf, 'paperposition', [ llx lly fpos(3) fpos(4) ] );
        end

        set(gcf,'paperunits',origPunits);   
        set(gcf,'units',origFunits);