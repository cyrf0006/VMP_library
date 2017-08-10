function [] = fig2gui(FigHndls,guiName)
%
% fig2gui.m--Given the handle to one or more figure windows containing
% uicontrol objects, fig2gui.m generates all the m-files required for a
% functioning GUI incorporating those figure windows.
%
% Input arguments are FigHndls, a vector of handles to the GUI figure
% windows, and guiName, a string giving the name of the GUI application to
% be created.
%
% Several m-files are created by fig2gui.m:
%
% a) The main GUI m-file--If guiName is 'mygui', the main m-file will have
% the name 'mygui.m'. The user invokes the GUI application, once it is
% created, by typing 'mygui'.
%
% b) "makeguiobjects" m-files--There will be one of these for each figure
% window in the GUI application.
%
% c) Uicontrol callback m-files--There will be one of these for each
% uicontrol object with a callback property in the GUI.
%
% d) The "getconstants m-file--This can be edited by the user to modify the
% GUI appearance (font size, figure window placement, etc.).
%
% e) The "start" m-file--Invoked when the GUI starts up.
%
% All m-files created by fig2gui.m will appear in the current directory.
%
% Before fig2gui.m can be run, it is necessary to create the GUI figure
% window with the desired uicontrol objects (pushbuttons, checkboxes,
% etc.). This is most easily done using Matlab's "guide" utility. 
%
% Each figure window to be included in the GUI must have a non-empty and
% unique 'Tag' property. Likewise, uicontrol objects should have unique
% 'Tag' properties (uicontrol objects with empty tags will still appear in
% the completed GUI, but no callback function will be created for them). 
%
% Only certain graphical objects and their properties are recognised by
% fig2gui.m. There is no point, for example, in using the guide utility
% to assign a figure window's 'PaperSize' property, since fig2gui.m will
% ignore it when generating its "makeguiobjects" code. The same goes for
% callback properties for uicontrols. Type 'fig2mfile' at the command line
% to see a list of properties supported by fig2gui.m.
%
% STEPS IN USING FIG2GUI.M:
% 
% 1) Create the desired uicontrol objects in one or more figure windows.
% The easiest way to do this is to use Matlab's "guide" utility.
%
% 2) In the desired directory, type "fig2gui(FigHndls,guiName)", where
% FigHndls is a vector of figure window handles and guiName is a string
% containing the name of the GUI you are creating. 
%
% 3) Edit the "start" m-file just created by fig2gui.m. If your GUI name
% was 'mygui', this file will have the name 'mygui_start.m'. In this file,
% you can make modifications to the GUI's appearance, changing object
% properties that are not recognised by fig2gui.m, or creating graphical
% objects that fig2gui.m ignores. For example, you might set the 'Box'
% property of a set of axes to 'on', and create some line objects the GUI
% needs. By default, all the GUI's figure windows are created whenever the
% GUI starts up. Override this by deleting the corresponding
% "makeguiobjects" lines in the start function. 
%
% 4) Edit the "getconstants" m-file just created by fig2gui.m. If your
% GUI name was 'mygui', this file will have the name
% 'mygui_getconstants.m'. Editing the values in this file will change
% figure scaling, font size, etc. This will generally not be necessary if
% you are using the GUI on the machine you designed it on, but is useful if
% you want to run it on a different computer or with a different screen
% resolution, etc.
%
% 5) Edit the callback m-files just created by fig2gui.m so that your
% uicontrols actually perform some useful function.
%
% 6) Close your figure windows and invoke the GUI. If your GUI name was
% 'mygui', invoke it by typing 'mygui' at the command line.
%
% Syntax: fig2gui(FigHndls,guiName)
%
% e.g.,   % Create figures with uicontrols.
%         f1 = figure('Tag','guifig01'); 
%         uicontrol('Tag','mypushbutton','Style','pushbutton','string','GO'); 
%         f2 = figure('Tag','guifig02'); 
%         uicontrol('Tag','myslider','Style','slider'); 
%     
%         % Generate GUI m-files.
%         fig2gui([f1 f2],'testgui');
% 
%         % Run the GUI.
%         close([f1 f2]);
%         testgui;

% Developed in Matlab 6.1.0.450 (R12.1) on Linux.
% Kevin Bartlett(kpb@hawaii.edu), 2003/06/13, 17:21
%------------------------------------------------------------------------------

if nargin ~= 2
   error([mfilename '.m--Incorrect number of input arguments.']);
end % if

% Check that all figures have unique, non-empty Tag properties.
FigHndls = FigHndls(:);
NumFigs = length(FigHndls);
TagCell = get(FigHndls,'Tag');

if isempty(TagCell)
   error([mfilename '.m--Figure ''Tag'' properties must not be empty.']);
end % if

% TagCell will be a character array if it only has one element.
if ~iscell(TagCell)
   TagData = TagCell;
   TagCell = {};
   TagCell{1} = TagData;
end % if

if size(strvcat(TagCell{:}),1)<NumFigs
   error([mfilename '.m--Figure ''Tag'' properties must not be empty.']);
end % if

if size(unique(TagCell),1)<NumFigs
   error([mfilename '.m--Figure ''Tag'' properties must be unique.']);
end % if

% Assemble filenames for output.

% ...Need names of uicontrol objects to make their callback functions.
AllCtrlList = findobj(FigHndls,'Type','uicontrol');

% ......Exclude uicontrol objects that do not normally have callbacks.
AllCtrlStyles = get(AllCtrlList,'Style');
ExcludeIndex = [strmatch('frame',AllCtrlStyles); strmatch('text',AllCtrlStyles)];
KeepIndex = setdiff([1:length(AllCtrlStyles)],ExcludeIndex);
CallbackCtrlList = AllCtrlList(KeepIndex);

if isempty(CallbackCtrlList)
   error([mfilename '.m--Cannot make GUI. There are no uicontrol objects in the figure window(s).']);
end % if

OutNames = fig2gui_filenames(guiName,CallbackCtrlList);

% See if any of the files exist already.
AllNames = cellstr(strvcat(OutNames.mainFileName,OutNames.getconstantsFileName,OutNames.startFileName,strvcat(OutNames.makeguiobjectsFileNames{:}),strvcat(OutNames.CallBackNames{:})));
FileExists = 0;

for FileCount = 1:length(AllNames)
   
   ThisFileName = fullfile(pwd,[AllNames{FileCount} '.m']);
   
   if exist(ThisFileName,'file')==2      
      FileExists = 1;
   end % if
   
end % for

if FileExists == 1
   
   disp(['Files for a GUI named ''' guiName ''' already exist.']);
   r = input(['Overwrite? (Y/N)   '],'s');
   
   if strcmp(lower(r),'y') ~= 1
      disp([mfilename ' execution halted. GUI m-files NOT created.']);
      return;
   end % if
   
end % if

% Create the "makeguiobjects" file for each figure.
ForceOverWrite = 1;

for FigCount = 1:NumFigs
   
   ThisFig = FigHndls(FigCount);
   
   ThisFigOutNames = fig2gui_filenames(guiName,ThisFig);
   ThisFileName = [ThisFigOutNames.makeguiobjectsFileNames{1} '.m'];
   
   % Create the "makeguiobjects" function.   
   fig2mfile(ThisFig,ThisFileName,ForceOverWrite);
   
end % for each figure.

% Create the "getconstants" function.
mk_getconstants_fcn(FigHndls,guiName);

% Create the "start" initialisation function.
mk_start_fcn(FigHndls,guiName);

% Create the main switchyard function.
mk_main_fcn(guiName);

% Create the callback function for each uicontrol object.
mk_callback_fcns(CallbackCtrlList,guiName);

% Copy files required by new GUI application to the current directory.
% N.B., there is an apparent bug that causes Matlab's platform-independent
% copyfile.m program to fail on my unix and linux boxes, so will use "unix"
% copy command where possible.
mFilesToCopy = {'parvalpairs' 'resizeguifig' 'repositiongui' 'resizeguifonts'};

CompStr = computer;

for mFileCount = 1:length(mFilesToCopy)

    ThisFileBaseName = mFilesToCopy{mFileCount};
    ThisFileFullName = which(ThisFileBaseName);
    ThisFileOutName = fullfile(pwd,[ThisFileBaseName '.m']);

    if ~isequal(ThisFileFullName,ThisFileOutName)

        success = 0;

        try
            [success,CopyMssg,CopyMssgID] = copyfile(ThisFileFullName,ThisFileOutName);
        catch
            success = 0;
        end % if

        % If copyfile.m didn't work, try unix copy commands, if available.
        if success == 0

            if strcmp(CompStr,'GLNX86')
                [status,result] = unix(['/bin/cp -p ' ThisFileFullName ' ' ThisFileOutName]);
            elseif strcmp(CompStr,'SOL2')
                [status,result] = unix(['/usr/bin/cp -p ' ThisFileFullName ' ' ThisFileOutName]);
            end % if

            success = ~status;
            
        end % if success == 0

        % If success is 0 for this file, then program fails.
        if success == 0
            error([mfilename '.m--Failed to copy all file(s) required for GUI to function.']);            
        end % if
        
    end % if ~isequal(ThisFileFullName,ThisFileOutName)

end % for

   

