function [subjects] = concatenate_take_files(Session)

dirArray = dir(fullfile(Session.spath, Session.maze.name));