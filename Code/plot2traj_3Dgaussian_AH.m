function plot2traj_3Dgaussian_AH(data, smooth_window, d_legend, d_marker_loc, d_marker_size, d_marker_legend, l_color, l_opacity, l_width, p_color, p_size, p_freq, a_xlim, a_ylim, a_xlabel, a_ylabel, a_title, a_fsize)
% Kanha Batra's neural trajectory script, adaptation to other smoothing
% Aniek
%    

% Purpose: plot 2 trajectories in the same figure using gaussian
    %          smoothening
    %
    % Inputs:
    % data : 1x6 cell array
    %     each element is a double array with dimensions 2xtimexcases 
    %     (denoted by T and C respectively)
    % d_smooth : double between 0 and 1
    %     smoothening factor for the data
    % d_legend : 1x6 cell array
    %     each element is a string for the legend entry of the
    %     corresponding trajectory
    % d_marker_loc : 1x3 double array
    %     each element represents the onset of an event in the time series,
    %     eg baseline, cue onset and trial end; each value must lie between
    %     0 and T
    % d_marker_size : double greater than 0
    %     size of the event marker
    % d_marker_legend : 1x3 cell array
    %     each element is a string for the name of the event
    % l_color : 1x6 cell array
    %     each element is a 1x3 double array for the color of the
    %     corresponding trajectory; each element of the array is between 0
    %     and 255, i.e., the RGB values
    % l_opacity : double between 0 and 1
    %     opacity for the line plot
    % l_width : double greater than 0
    %     width of the line plot
    % p_color : 1x6 cell array
    %     each element is a 1x3 double array for the color of the
    %     corresponding time bin markers; each element of the array is
    %     between 0 and 255, i.e., the RGB values
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
    
    % Code:
    d1 = smoothdata(mean(data{1}, 3), 2, 'gaussian', smooth_window);
    d2 = smoothdata(mean(data{2}, 3), 2, 'gaussian', smooth_window);
   
    %figure();
    hold on;
    p1 = plot3(d1(1, :), d1(2, :), d1(3,:), 'DisplayName', d_legend{1});
    p1.Color(1: 3) = l_color{1}; p1.Color(4) = l_opacity; p1.LineWidth = l_width;
    p1.Marker = '.'; p1.MarkerFaceColor = p_color{1}; p1.MarkerEdgeColor = p_color{1}; p1.MarkerIndices = [1: p_freq: size(d1, 2)]; p1.MarkerSize = p_size;
    e11 = scatter3(d1(1, d_marker_loc(1)), d1(2, d_marker_loc(1)), d1(3, d_marker_loc(1)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e12 = scatter3(d1(1, d_marker_loc(2)), d1(2, d_marker_loc(2)), d1(3, d_marker_loc(2)), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    e13 = scatter3(d1(1, d_marker_loc(3)), d1(2, d_marker_loc(3)), d1(3, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    
    p2 = plot3(d2(1, :), d2(2, :),d2(3, :), 'DisplayName', d_legend{2});
    p2.Color(1: 3) = l_color{2}; p2.Color(4) = l_opacity; p2.LineWidth = l_width;
    p2.Marker = '.'; p2.MarkerFaceColor = p_color{2}; p2.MarkerEdgeColor = p_color{2}; p2.MarkerIndices = [1: p_freq: size(d2, 2)]; p2.MarkerSize = p_size;
    e21 = scatter3(d2(1, d_marker_loc(1)), d2(2, d_marker_loc(1)), d2(3, d_marker_loc(1)), d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'HandleVisibility', 'off');
    e22 = scatter3(d2(1, d_marker_loc(2)), d2(2, d_marker_loc(2)), d2(3, d_marker_loc(2)), d_marker_size, [0, 0, 0]/255, '>', 'filled', 'HandleVisibility', 'off');
    e23 = scatter3(d2(1, d_marker_loc(3)), d2(2, d_marker_loc(3)), d2(3, d_marker_loc(3)), d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'HandleVisibility', 'off');
    
    e1 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 0]/255, 'o', 'filled', 'Displayname', d_marker_legend{1});
    e2 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 0]/255, '>', 'filled', 'Displayname', d_marker_legend{2});
    e3 = scatter3(NaN, NaN, NaN, d_marker_size, [0, 0, 0]/255, 'square', 'filled', 'Displayname', d_marker_legend{3});
    
    hold off;
    %legend();
    if nargin > 12
        xlim(a_xlim);
        ylim(a_ylim);
        xlabel(a_xlabel);
        ylabel(a_ylabel);
        title(a_title);
        set(findall(gcf,'-property','FontSize'),'FontSize',a_fsize)
    end
end
