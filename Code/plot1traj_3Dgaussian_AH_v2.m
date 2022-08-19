function plot1traj_3Dgaussian_AH_v2(data, drinking_data, smooth_window, d_legend, d_marker_size, d_marker_legend, l_color, l_opacity, l_width, p_color, p_size, p_freq, a_xlim, a_ylim, a_xlabel, a_ylabel, a_title, a_fsize)
% Kanha Batra's neural trajectory script, adaptation to full hour session trajectory colored by events by Aniek 
%
    % Purpose: plot 1 trajectories in the same figure using gaussian
    %          smoothening
    %
    % Inputs:
    % data : 1x1 cell array
    %     each element is a double array with dimensions 2xtimexcases 
    %     (denoted by T and C respectively)
    % drinking_data: Binary information on drinking (4xtimepoints matrix)
    % smooth_window : Array of 2 numbers, frames used for smoothing
    %     Always put 0 as second to not use future smoothing
    % d_legend : 1x1 cell array
    %     each element is a string for the legend entry of the
    %     corresponding trajectory
    % d_marker_size : double greater than 0
    %     size of the event marker
    % d_marker_legend : 1x3 cell array
    %     each element is a string for the name of the event
    % l_color : 1x1 cell array
    %     each element is a 1x3 double array for the color of the
    %     corresponding trajectory; each element of the array is between 0
    %     and 255
    % l_opacity : double between 0 and 1
    %     opacity for the line plot
    % l_width : double greater than 0
    %     width of the line plot
    % p_color : 1x1 cell array
    %     each element is a 1x3 double array for the color of the
    %     corresponding time bin markers; each element of the array is
    %     between 0 and 255  
    % p_size : double greater than 0
    %     size of time bin markers
    % p_freq : double between 0 and T
    %     frequency of markers, eg if p_freq = 5, every 5th bin will be
    %     marked
    %
    % Optional Input:
    % a_xlim : 1x2 double array
    %     x-axis limits
    % a_ylim : 1x2 double array
    %     y-axis limits
    % a_xlabel : string
    %     title of the x-axis units
    % a_ylabel : string
    %     title of the y-axis units
    % a_title : string
    %     title of the overall plot
    % a_fsize : double
    %     font size for the axes titles, axes labels, figure title and
    %     legend
    %
    % Output:
    % <none> : a new figure is constructed
    
    % Last updated: 08/19/2022 Aniek van Hoek
    
    % Code:
    d1 = smoothdata(data{1}, 2, 'gaussian', smooth_window);
   
    hold on;
    p1 = plot3(d1(1, :), d1(2, :), d1(3, :), 'DisplayName', d_legend{1});
    p1.Color(1:3) = l_color{1}; p1.Color(4) = l_opacity; p1.LineWidth = l_width;
    p1.Marker = '.'; p1.MarkerFaceColor = p_color{1}; p1.MarkerEdgeColor = p_color{1}; p1.MarkerIndices = [1: p_freq: size(d1, 2)]; p1.MarkerSize = p_size;
    
    start = scatter3(d1(1, 1), d1(2, 1), d1(3, 1), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    finish = scatter3(d1(1, end), d1(2, end), d1(3, end), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
         
    start2 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 0]/255, '>', 'filled', 'Displayname', d_marker_legend{1});
	finish2 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'Displayname', d_marker_legend{2});
    
    
    % Ethanol bouts
    e11 = scatter3(d1(1, logical(drinking_data(4,:))), d1(2, logical(drinking_data(4,:))), d1(3, logical(drinking_data(4,:))), d_marker_size, [255, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e1 = scatter3(NaN, NaN, NaN, d_marker_size, [255, 0, 0]/255, 'o', 'filled', 'Displayname', d_marker_legend{3});
    
    % Water bouts
    e12 = scatter3(d1(1, logical(drinking_data(3,:))), d1(2, logical(drinking_data(3,:))), d1(3, logical(drinking_data(3,:))), d_marker_size, [0, 0, 255]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e2 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 255]/255, 'o', 'filled', 'Displayname', d_marker_legend{4});
    

    % Cue indices
    e13 = scatter3(d1(1, logical(drinking_data(2,:))), d1(2, logical(drinking_data(2,:))), d1(3, logical(drinking_data(2,:))), d_marker_size, [0, 128, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e3 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 128, 0]/255, 'o', 'filled', 'Displayname', d_marker_legend{5});
        
    

    hold off;
    if nargin > 12
        xlim(a_xlim);
        ylim(a_ylim);
        xlabel(a_xlabel);
        ylabel(a_ylabel);
        title(a_title);
        set(findall(gcf,'-property','FontSize'),'FontSize',a_fsize)
    end
end
