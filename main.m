clc;
close all;
clear all;  % Clear cached classes
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.

% Setup environment
if count(py.sys.path,pwd) == 0
    insert(py.sys.path,int32(0),pwd);
end
setenv('TCL_LIBRARY', 'C:/Python312/tcl/tcl8.6');
setenv('TK_LIBRARY', 'C:/Python312/tcl/tk8.6');

import modules.*;

% workspace;  % Make sure the workspace panel is showing.
% start parallel pool
% parpool(4);
% TODO: Time stamp
% TODO: Provision to combine stacks
% TODO: Time stamp from OCR
% TODO: Plot all timestamps

searchWindow = 50;
function_list = {'Change Drive', 'Plot all Gr', 'Plot all TimeStamps'};
get_time = [];path = [];
is_playing = false;
speed = 1;
logs = {}; % Initialize logs list
stack_info = struct();
skip_alignment = false;
forced = false;

app = app(stack_paths);
goto_callback();
