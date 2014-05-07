classdef EMApexFloat2
    % This is a container class for all data relating to a single EM-APEX
    % float including all profiles. It will also (eventually include
    % methods for plotting using multiple profile data).
    
    properties
        ID          % Float name.
        pfls        % Array of Profile classes.
        nPfls       % Number of profiles.
        hpids       % Array of hpid numbers.
        MLT_refs    % Array of profile start times.
        lats        % Array of latitudes for corresponding hpids.
        lons        % Array of longitudes for corresponding hpids.
        
        % Calculated.
        LATS        % Array of interpolated latitudes for corresponding hpids.
        LONS        % Array of interpolated longtidues for corresponding hpids.
        DIST        % Array of distances from first profile calculated using LATS and LONS.
        U           % Array of eastward velocities calculated using LATS and LONS.
        V           % Array of northward velocities calculated using LATS and LONS.
        
        % Taken from the allprof##.mat file. See allprofs_report.pdf for
        % more information.
        allProfFNs = {'u_gps','u_sfc','ubs','v_gps','v_sfc','vbs'};
        u_gps
        u_sfc
        ubs
        v_gps
        v_sfc
        vbs
        
    end % properties
    
    methods
        function obj = EMApexFloat2(ID, path_to_file)

            disp('grrr')
            
        end % function 
        
        function [pflsSubset, validHpidIndx, validHpidIndxBool] = ...
                getPflSubset(obj, hpidIndx)
            % USAGE:
            %  getPflSubset(hpidIndx)
            %  
            % DESCRIPTION:
            %  Return profile(s) with given hpid numbers.
            %
            % INPUT:
            %  hpidIndx = array of half profile ID numbers.
            %
            % OUTPUT:
            %  pflsSubset = Array of Profiles. If input is a single
            %  number, output will be a single Profile instance. If input
            %  is a single number and Profile hpid does not exist, output
            %  will will be an empty Profile, beware. 
            %  validHpidIndx = Array containing the hpid numbers of the
            %  floats that were successfully found. 
            %  validHpidIndxBool = Boolean array of valid hpid indices.
            %
            % EXAMPLE:
            %  % Retrieve Profile 30.
            %  EMF.getPflSubset(30)
            %  % Retrieve Profiles 10, 12, 14 ... 200.
            %  EMF.getPflSubset(10:2:200)
            
            validHpidIndxBool = ismember(obj.hpids, hpidIndx);
            validHpidIndx = obj.hpids(validHpidIndxBool);
            pflsSubset = obj.pfls(validHpidIndxBool);
            
        end % function
        
        function [h, ax1, ax2] = section(obj, hpidIndx, xVar, zVar, ...
                griddedVar, plotfunc, anom)
            % USAGE:
            %  section(hpidIndx, xVar, zVar, griddedVar, plotfunc, anom)
            %  
            % DESCRIPTION:
            %  Plot a property measured by the float on a section. The z
            %  axis variable zVar can be specified, for example, 'z' is
            %  depth. The horizontal axis has to be specified as one of the
            %  variables in the EMApexFloat class e.g. LON for longitude,
            %  DIST for distance, MLT_refs for time. The scattered data is
            %  gridded against these two axes using the griddata function.
            %
            % INPUT:
            %  hpidIndx = array of half profile ID numbers to plot.
            %  xVar = property name of the xVar, e.g. 'DIST'.
            %  zVar = property name of the zVar. e.g. 'z' or 'P'.
            %  griddedVar = property name of the gridded variable, e.g.
            %  'CT'.
            %  plotfunc = handle to a suitable 2-D plotting function that
            %  can handle grid inputs e.g. @pcolor or @contourf.
            %  anom = if true will subtract row averaged data before
            %  plotting, resulting in an anomaly plot.
            %
            % OUTPUT:
            %  h = figure handle.
            %  ax1 = first axis handle (where the data is plotted)
            %  ax2 = second axis handle (where the profile points are
            %  denoted).

            [pflsSubset, hpidIndxs, indxBools] = obj.getPflSubset(hpidIndx);
            N = length(pflsSubset);
            
            z = [];
            x = [];
            prop = [];
            
            % X axis variable for each chosen half profile.
            X = obj.(xVar)(indxBools);
            
            % The following conditions check whether the data is measured
            % at the time of the efp or the cdt, and interpolates zVar
            % to the time of measurements if necessary. 
            testPfl = pflsSubset(1);
            if length(testPfl.(griddedVar)) == length(testPfl.(zVar))
                
                i = 1;
                for pfl = pflsSubset'

                    z = [z; pfl.(zVar)];
                    x = [x; X(i)*ones(size(pfl.(zVar)))];
                    prop = [prop; pfl.(griddedVar)];
                    i = i + 1;
                    
                end % for
                
            elseif length(testPfl.(griddedVar)) < length(testPfl.(zVar)) && ...
                    length(testPfl.(griddedVar)) == length(testPfl.efp_mlt) && ...
                    length(testPfl.(zVar)) == length(testPfl.ctd_mlt)
                
                i = 1;
                for pfl = pflsSubset'

                    z = [z; ...
                        interp1(pfl.ctd_mlt, pfl.(zVar), pfl.efp_mlt, ...
                        'linear', 'extrap')];
                    x = [x; X(i)*ones(size(pfl.(griddedVar)))];
                    prop = [prop; pfl.(griddedVar)];
                    i = i + 1;
                    
                end % for
                
            else
                error('Have you chosen an appropriate variables to section?')
            end % if
            
            % 
            [Xg, Pg] = meshgrid(linspace(min(x), max(x), 5*N), ...
                linspace(min(z), max(z), 500));
            prop_grid = griddata(x, z, prop, Xg, Pg);
             
            if anom
                anom_grid = prop_grid - ...
                    repmat(nanmean(prop_grid, 2), 1, size(prop_grid, 2));
                h = plotfunc(Xg, Pg, anom_grid);
            else
                h = plotfunc(Xg, Pg, prop_grid);
            end

            % Add some useful but basic information. In the interest of
            % simplicity, more complex figure attributes should be added
            % later (e.g. colour bars).
            xlabel(xVar);
            ylabel(zVar);
            title([obj.ID ' ' griddedVar]);
            ax1 = gca();
            set(ax1, 'Box', 'off');
            ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none');
            set(ax2, 'XAxisLocation', 'top');
            set(ax2, 'XLim', get(ax1, 'XLim'));
            set(ax2, 'XTick', X, 'YTick', []);
            set(ax2, 'XTickLabel', hpidIndxs, 'YTickLabel', []);
            set(ax2, 'YAxisLocation', 'right')
            set(ax2, 'YLim', get(ax1, 'YLim'))
            set(ax2, 'YTick', get(ax1, 'Ytick')) 

        end % function
        
        function [h] = plot_track(obj, hpidIndx, annotate)
            % USAGE:
            %  plot_track(hpidIndx)
            %  
            % DESCRIPTION:
            %  Plot the trajectory of the float for the given profiles. It
            %  uses GPS positions determined by interpolation from the GPS
            %  file (more trustworthy) and also that derived from the 'vel'
            %  files (less trustworthy). 
            %
            %  Note: I have commented out the less trustworthy positions. 
            %
            % INPUT:
            %  hpidIndx = array of half profile ID numbers to grid.
            %  annotate = true/false flag. 
            %
            % OUTPUT:
            %  None
            
            [~, hpidIndxs, indxBool] = obj.getPflSubset(hpidIndx);
            
%             lats = obj.lats(indxBool); lons = obj.lons(indxBool);
            LATS = obj.LATS(indxBool); LONS = obj.LONS(indxBool);
            
            % Extract grid...
            margin = 1;
            latlim = [min(LATS)-margin, max(LATS)+margin];
            lonlim = [min(LONS)-margin, max(LONS)+margin];
            
            [ele, elelat, elelon] = mygrid_sand([latlim lonlim], 1);
            ele(ele > 0) = NaN;
            
            h = figure;
            hold on
            m_proj('lambert','lon',lonlim,'lat',latlim); 
            m_grid('box','on','color','k','linewidth',1,'fontsize',12);
            m_pcolor(elelon - 360, elelat, ele); 
            colormap Jet;
            shading flat; 
            c = colorbar;
            ylabel(c, 'Depth (m)')
%             m_plot(lons, lats, 'k.');
            m_plot(LONS, LATS, 'ro');
            
            if annotate
            % Annotate profile with a number.
                for i = 1:length(LATS)
%                     m_text(lons(i), lats(i), num2str(hpidIndxs(i)))
                    m_text(LONS(i), LATS(i), num2str(hpidIndxs(i)), 'Color', 'r')
                end
            end

            hold off
            
        end % function
        
        function [xg zg propg] = gridVar(obj, hpidIndx, zVar, zLvls, ...
                griddedVar)
            % USAGE:
            %  gridVar(hpidIndx, zVar, zLvls, griddedVar)
            %  
            % DESCRIPTION:
            %  Produce a grid of some property at the given levels of some
            %  other property (e.g. depth) and for the specified half
            %  profiles. 
            %
            % INPUT:
            %  hpidIndx = array of half profile ID numbers to grid.
            %  zVar = property name of the zVar. e.g. 'z' or 'P'.
            %  griddedVar = property name of the gridded variable, e.g.
            %  'CT' or 'pdens'.
            %
            % OUTPUT:
            %  xg = grid of arbitrary counter values whos maximum is equal
            %  to the total number of profiles included in the grid.
            %  zg = grid of the zVar.
            %  propg = grid of the property specified with griddedVar.
            
            pflsSubset = obj.getPflSubset(hpidIndx);
            
            i = 0;
            for pfl = pflsSubset'
                i = i + 1;
                propg(:,i) = pfl.interp_var(griddedVar, zVar, zLvls);
            end % for
            
            [xg, zg] = meshgrid(1:i, zLvls);
            
        end % function
        
        function [] = calc_vert_vel(obj, up_params, down_params)
            % Calls the calc_vert_vel function for each profile. Assumes
            % that separate parameters are required for up and down
            % profiles.
            
            even_pfls = obj.getPflSubset(obj.hpids(rem(obj.hpids, 2) == 0));
            odd_pfls = obj.getPflSubset(obj.hpids(rem(obj.hpids, 2) == 1));
            
            for pfl = even_pfls'
                
                pfl.calc_vert_vel(up_params)
                
            end % for
            
            for pfl = odd_pfls'
                
                pfl.calc_vert_vel(down_params)
                
            end % for
            
        end % function
        
        function [] = plot_subset(obj, func, hpids, propx, propy)
            
            [pfls, validHpids, ~] = obj.getPflSubset(hpids);
            N = length(pfls);
            
            figure;
            hold all
            
            colours = [linspace(0,1,N)' zeros(N,1) linspace(1,0,N)'];
            
            i = 1;
            for pfl = pfls'
                
                pfl.plot(func, propx, propy, 'color', colours(i, :))
                i = i + 1;
                
                
            end % for
            
            legend(num2str(validHpids))
            hold off
            
        end % function
        
    end % methods
    
end % classdef