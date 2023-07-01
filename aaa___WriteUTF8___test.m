function pass_part_fail=aaa___WriteUTF8___test(varargin)
%
%Pass:    passes all tests
%Partial: [no partial passing condition]
%Fail:    fails any test
pass_part_fail='pass';

if nargin==0,RunTestHeadless = false;else,RunTestHeadless = true;end

% Create a target folder name.
checkpoint('aaa___WriteUTF8___test','mexname')
[ignore,base] = mexname('WriteUTF8_test'); %#ok<ASGLU>
[p_,f_] = fileparts(tempname);
p = [p_ filesep base f_];

% Generate some rich text to be written.
checkpoint('aaa___WriteUTF8___test','ch_r')
str = ch_r([100 169 233 239 8252]);
str_uint32 = uint32([100 169 233 239 8252]);

ErrorFlag = false;
try ME=[];
    checkpoint('aaa___WriteUTF8___test','makedir')
    makedir(p); % Create the folder.
    WriteUTF8([p filesep 'with_BOM.txt']   ,str,'BOM',true )
    WriteUTF8([p filesep 'without_BOM.txt'],str,'BOM',false)
catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
    ErrorFlag = true;
end

if ~ErrorFlag
    try
        % These three files should all have the same contents.
        % x: Write baseline, but with UTF-32.
        % y: Write with a text mode permission.
        % z: Attempt to write a BOM on the second pass.
        WriteUTF8([p filesep 'x.txt'],str_uint32,'BOM',false,'permission','w')
        WriteUTF8([p filesep 'x.txt'],str_uint32,'BOM',false,'permission','a')
        WriteUTF8([p filesep 'y.txt'],str       ,'BOM',false,'permission','wt')
        WriteUTF8([p filesep 'y.txt'],str       ,'BOM',false,'permission','at')
        WriteUTF8([p filesep 'z.txt'],str       ,'BOM',false,'permission','w')
        WriteUTF8([p filesep 'z.txt'],str       ,'BOM',true ,'permission','a')
        contents = {...
            WriteUTF8_helper___read_binary([p filesep 'x.txt']) ,...
            WriteUTF8_helper___read_binary([p filesep 'y.txt']) ,...
            WriteUTF8_helper___read_binary([p filesep 'z.txt']) };
        if ~isequal(contents{1},contents{2})
            error('permission test failed: files should be equal.')
        end
        if ~isequal(contents{1},contents{3})
            error('BOM test failed: files should be equal.')
        end
        % Test succeeded, so clean up
        delete([p filesep 'x.txt']);
        delete([p filesep 'y.txt']);
        delete([p filesep 'z.txt']);
    catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
        ErrorFlag = true;
    end
end

checkpoint('aaa___WriteUTF8___test','ifversion')
persistent isOctave,if isempty(isOctave),isOctave=ifversion('<',0,'Octave','>',0);end
if RunTestHeadless
    % During normal operation (if in headless mode) delete the folder and the files.
    try
        if isOctave,oldval=confirm_recursive_rmdir(false);end
        rmdir(p,'s');
        if isOctave,confirm_recursive_rmdir(oldval);end
    catch
    end
else
    % Open the directory so the user can check the files.
    checkpoint('aaa___WriteUTF8___test','OpenFileExplorer')
    OpenFileExplorer(p)
end

SelfTestFailMessage = '';
% Run the self-validator function(s).
SelfTestFailMessage = [SelfTestFailMessage SelfTest__PatternReplace];
SelfTestFailMessage = [SelfTestFailMessage SelfTest__error_];
checkpoint('read');
if ~isempty(SelfTestFailMessage) || ErrorFlag
    if nargout==1
        pass_part_fail='fail';
    else
        if ~isempty(SelfTestFailMessage)
            error('Self-validator functions returned these error(s):\n%s',SelfTestFailMessage)
        else
            error('Test failed with message:\n    %s',ME.message)
        end
    end
end
disp(['tester function ' mfilename ' finished '])
if nargout==0,clear,end
end
function bin=WriteUTF8_helper___read_binary(filename)
% Read the file to a byte stream.
fid = fopen(filename,'r');
bin = fread(fid,inf,'uint8=>uint8');
fclose(fid);
end
function out=bsxfun_plus(in1,in2)
%Implicit expansion for plus(), but without any input validation.
persistent type
if isempty(type)
    checkpoint('bsxfun_plus','hasFeature')
    type = ...
        double(hasFeature('ImplicitExpansion')) + ...
        double(hasFeature('bsxfun'));
end
if type==2
    % Implicit expansion is available.
    out = in1+in2;
elseif type==1
    % Implicit expansion is only available with bsxfun.
    out = bsxfun(@plus,in1,in2);
else
    % No implicit expansion, expand explicitly.
    sz1 = size(in1);
    sz2 = size(in2);
    if min([sz1 sz2])==0
        % Construct an empty array of the correct size.
        sz1(sz1==0) = inf;sz2(sz2==0) = inf;
        sz = max(sz1,sz2);
        sz(isinf(sz)) = 0;
        % Create an array and cast it to the correct type.
        out = feval(str2func(class(in1)),zeros(sz));
        return
    end
    in1 = repmat(in1,max(1,sz2./sz1));
    in2 = repmat(in2,max(1,sz1./sz2));
    out = in1+in2;
end
end
function out=ch_r(in)
checkpoint('ch_r','unicode_to_char','CharIsUTF8')
out = unicode_to_char(in,~CharIsUTF8);
end

function c=char2cellstr(str,LineEnding)
% Split char or uint32 vector to cell (1 cell element per line). Default splits are for CRLF/CR/LF.
% The input data type is preserved.
%
% Since the largest valid Unicode codepoint is 0x10FFFF (i.e. 21 bits), all values will fit in an
% int32 as well. This is used internally to deal with different newline conventions.
%
% The second input is a cellstr containing patterns that will be considered as newline encodings.
% This will not be checked for any overlap and will be processed sequentially.

returnChar = isa(str,'char');
str = int32(str); % Convert to signed, this should not crop any valid Unicode codepoints.

if nargin<2
    % Replace CRLF, CR, and LF with -10 (in that order). That makes sure that all valid encodings
    % of newlines are replaced with the same value. This should even handle most cases of files
    % that mix the different styles, even though such mixing should never occur in a properly
    % encoded file. This considers LFCR as two line endings.
    if any(str==13)
        checkpoint('char2cellstr','PatternReplace')
        str = PatternReplace(str,int32([13 10]),int32(-10));
        str(str==13) = -10;
    end
    str(str==10) = -10;
else
    for n=1:numel(LineEnding)
        checkpoint('char2cellstr','PatternReplace')
        str = PatternReplace(str,int32(LineEnding{n}),int32(-10));
    end
end

% Split over newlines.
newlineidx = [0 find(str==-10) numel(str)+1];
c=cell(numel(newlineidx)-1,1);
for n=1:numel(c)
    s1 = (newlineidx(n  )+1);
    s2 = (newlineidx(n+1)-1);
    c{n} = str(s1:s2);
end

% Return to the original data type.
if returnChar
    for n=1:numel(c),c{n} =   char(c{n});end
else
    for n=1:numel(c),c{n} = uint32(c{n});end
end
end
function tf=CharIsUTF8
% This provides a single place to determine if the runtime uses UTF-8 or UTF-16 to encode chars.
% The advantage is that there is only 1 function that needs to change if and when Octave switches
% to UTF-16. This is unlikely, but not impossible.
persistent persistent_tf
if isempty(persistent_tf)
    checkpoint('CharIsUTF8','ifversion')
    if ifversion('<',0,'Octave','>',0)
        % Test if Octave has switched to UTF-16 by looking if the Euro symbol is losslessly encoded
        % with char.
        % Because we will immediately reset it, setting the state for all warnings to off is fine.
        w = struct('w',warning('off','all'));[w.msg,w.ID] = lastwarn;
        persistent_tf = ~isequal(8364,double(char(8364)));
        warning(w.w);lastwarn(w.msg,w.ID); % Reset warning state.
    else
        persistent_tf = false;
    end
end
tf = persistent_tf;
end
function error_(options,varargin)
%Print an error to the command window, a file and/or the String property of an object.
% The error will first be written to the file and object before being actually thrown.
%
% Apart from controlling the way an error is written, you can also run a specific function. The
% 'fcn' field of the options must be a struct (scalar or array) with two fields: 'h' with a
% function handle, and 'data' with arbitrary data passed as third input. These functions will be
% run with 'error' as first input. The second input is a struct with identifier, message, and stack
% as fields. This function will be run with feval (meaning the function handles can be replaced
% with inline functions or anonymous functions).
%
% The intention is to allow replacement of every error(___) call with error_(options,___).
%
% NB: the function trace that is written to a file or object may differ from the trace displayed by
% calling the builtin error/warning functions (especially when evaluating code sections). The
% calling code will not be included in the constructed trace.
%
% There are two ways to specify the input options. The shorthand struct described below can be used
% for fast repeated calls, while the input described below allows an input that is easier to read.
% Shorthand struct:
%  options.boolean.IsValidated: if true, validation is skipped
%  options.params:              optional parameters for error_ and warning_, as explained below
%  options.boolean.con:         only relevant for warning_, ignored
%  options.fid:                 file identifier for fprintf (array input will be indexed)
%  options.boolean.fid:         if true print error to file
%  options.obj:                 handle to object with String property (array input will be indexed)
%  options.boolean.obj:         if true print error to object (options.obj)
%  options.fcn                  struct (array input will be indexed)
%  options.fcn.h:               handle of function to be run
%  options.fcn.data:            data passed as third input to function to be run (optional)
%  options.boolean.fnc:         if true the function(s) will be run
%
% Full input description:
%   print_to_con:
%      NB: An attempt is made to use this parameter for warnings or errors during input parsing.
%      A logical that controls whether warnings and other output will be printed to the command
%      window. Errors can't be turned off. [default=true;]
%      Specifying print_to_fid, print_to_obj, or print_to_fcn will change the default to false,
%      unless parsing of any of the other exception redirection options results in an error.
%   print_to_fid:
%      NB: An attempt is made to use this parameter for warnings or errors during input parsing.
%      The file identifier where console output will be printed. Errors and warnings will be
%      printed including the call stack. You can provide the fid for the command window (fid=1) to
%      print warnings as text. Errors will be printed to the specified file before being actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_obj:
%      NB: An attempt is made to use this parameter for warnings or errors during input parsing.
%      The handle to an object with a String property, e.g. an edit field in a GUI where console
%      output will be printed. Messages with newline characters (ignoring trailing newlines) will
%      be returned as a cell array. This includes warnings and errors, which will be printed
%      without the call stack. Errors will be written to the object before the error is actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_fcn:
%      NB: An attempt is made to use this parameter for warnings or errors during input parsing.
%      A struct with a function handle, anonymous function or inline function in the 'h' field and
%      optionally additional data in the 'data' field. The function should accept three inputs: a
%      char array (either 'warning' or 'error'), a struct with the message, id, and stack, and the
%      optional additional data. The function(s) will be run before the error is actually thrown.
%      [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_params:
%      NB: An attempt is made to use this parameter for warnings or errors during input parsing.
%      This struct contains the optional parameters for the error_ and warning_ functions.
%      Each field can also be specified as ['print_to_option_' parameter_name]. This can be used to
%      avoid nested struct definitions.
%      ShowTraceInMessage:
%        [default=false] Show the function trace in the message section. Unlike the normal results
%        of rethrow/warning, this will not result in clickable links.
%      WipeTraceForBuiltin:
%        [default=false] Wipe the trace so the rethrow/warning only shows the error/warning message
%        itself. Note that the wiped trace contains the calling line of code (along with the
%        function name and line number), while the generated trace does not.
%
% Syntax:
%   error_(options,msg)
%   error_(options,msg,A1,...,An)
%   error_(options,id,msg)
%   error_(options,id,msg,A1,...,An)
%   error_(options,ME)               %equivalent to rethrow(ME)
%
% Examples options struct:
%   % Write to a log file:
%   opts = struct;opts.fid = fopen('log.txt','wt');
%   % Display to a status window and bypass the command window:
%   opts = struct;opts.boolean.con = false;opts.obj = uicontrol_object_handle;
%   % Write to 2 log files:
%   opts = struct;opts.fid = [fopen('log2.txt','wt') fopen('log.txt','wt')];

persistent this_fun
if isempty(this_fun),this_fun = func2str(@error_);end

% Parse options struct, allowing an empty input to revert to default.
if isempty(options),options = struct;end
checkpoint('error_','parse_warning_error_redirect_options')
options                    = parse_warning_error_redirect_options(  options  );
checkpoint('error_','parse_warning_error_redirect_inputs')
[id,msg,stack,trace,no_op] = parse_warning_error_redirect_inputs( varargin{:});
if no_op,return,end
forced_trace = trace;
if options.params.ShowTraceInMessage
    msg = sprintf('%s\n%s',msg,trace);
end
ME = struct('identifier',id,'message',msg,'stack',stack);
if options.params.WipeTraceForBuiltin
    ME.stack = stack('name','','file','','line',[]);
end

% Print to object.
if options.boolean.obj
    msg_ = msg;while msg_(end)==10,msg_(end) = '';end % Crop trailing newline.
    if any(msg_==10)  % Parse to cellstr and prepend 'Error: '.
        checkpoint('error_','char2cellstr')
        msg_ = char2cellstr(['Error: ' msg_]);
    else              % Only prepend 'Error: '.
        msg_ = ['Error: ' msg_];
    end
    for OBJ=reshape(options.obj,1,[])
        try set(OBJ,'String',msg_);catch,end
    end
end

% Print to file.
if options.boolean.fid
    T = datestr(now,31); %#ok<DATST,TNOW1> Print the time of the error to the log as well.
    for FID=reshape(options.fid,1,[])
        try fprintf(FID,'[%s] Error: %s\n%s',T,msg,trace);catch,end
    end
end

% Run function.
if options.boolean.fcn
    if ismember(this_fun,{stack.name})
        % To prevent an infinite loop, trigger an error.
        error('prevent recursion')
    end
    ME_ = ME;ME_.trace = forced_trace;
    for FCN=reshape(options.fcn,1,[])
        if isfield(FCN,'data')
            try feval(FCN.h,'error',ME_,FCN.data);catch,end
        else
            try feval(FCN.h,'error',ME_);catch,end
        end
    end
end

% Actually throw the error.
rethrow(ME)
end
function [valid,filename]=filename_is_valid(filename)
% Check if the file name and path are valid (non-empty char or scalar string).
valid=true;
persistent forbidden_names
if isempty(forbidden_names)
    forbidden_names = {'CON','PRN','AUX','NUL','COM1','COM2','COM3','COM4','COM5','COM6','COM7',...
        'COM8','COM9','LPT1','LPT2','LPT3','LPT4','LPT5','LPT6','LPT7','LPT8','LPT9'};
end
if isa(filename,'string') && numel(filename)==1
    % Convert a scalar string to a char array.
    filename = char(filename);
end
if ~isa(filename,'char') || numel(filename)==0
    valid = false;return
else
    % File name is indeed a char. Do a check if there are characters that can't exist in a normal
    % file name. The method used here is not fool-proof, but should cover most use cases and
    % operating systems.
    [fullpath,fn,ext] = fileparts(filename); %#ok<ASGLU>
    fn = [fn,ext];
    if      any(ismember([char(0:31) '<>:"/\|?*'],fn)) || ...
            any(ismember(forbidden_names,upper(fn))) || ... % (ismember is case sensitive)
            any(fn(end)=='. ')
        valid = false;return
    end
end
end
function [str,stack]=get_trace(skip_layers,stack)
if nargin==0,skip_layers = 1;end
if nargin<2, stack = dbstack;end
stack(1:skip_layers) = [];

% Parse the ML6.5 style of dbstack (the name field includes full file location).
if ~isfield(stack,'file')
    for n=1:numel(stack)
        tmp = stack(n).name;
        if strcmp(tmp(end),')')
            % Internal function.
            ind = strfind(tmp,'(');
            name = tmp( (ind(end)+1):(end-1) );
            file = tmp(1:(ind(end)-2));
        else
            file = tmp;
            [ignore,name] = fileparts(tmp); %#ok<ASGLU>
        end
        [ignore,stack(n).file] = fileparts(file); %#ok<ASGLU>
        stack(n).name = name;
    end
end

% Parse Octave style of dbstack (the file field includes full file location).
checkpoint('get_trace','ifversion')
persistent isOctave,if isempty(isOctave),isOctave=ifversion('<',0,'Octave','>',0);end
if isOctave
    for n=1:numel(stack)
        [ignore,stack(n).file] = fileparts(stack(n).file); %#ok<ASGLU>
    end
end

% Create the char array with a (potentially) modified stack.
s = stack;
c1 = '>';
str = cell(1,numel(s)-1);
for n=1:numel(s)
    [ignore_path,s(n).file,ignore_ext] = fileparts(s(n).file); %#ok<ASGLU>
    if n==numel(s),s(n).file = '';end
    if strcmp(s(n).file,s(n).name),s(n).file = '';end
    if ~isempty(s(n).file),s(n).file = [s(n).file '>'];end
    str{n} = sprintf('%c In %s%s (line %d)\n',c1,s(n).file,s(n).name,s(n).line);
    c1 = ' ';
end
str = horzcat(str{:});
end
function tf=hasFeature(feature)
% Provide a single point to encode whether specific features are available.
persistent FeatureList
if isempty(FeatureList)
    checkpoint('hasFeature','ifversion')
    FeatureList = struct(...
        'HG2'              ,ifversion('>=','R2014b','Octave','<' ,0),...
        'ImplicitExpansion',ifversion('>=','R2016b','Octave','>' ,0),...
        'bsxfun'           ,ifversion('>=','R2007a','Octave','>' ,0),...
        'IntegerArithmetic',ifversion('>=','R2010b','Octave','>' ,0),...
        'String'           ,ifversion('>=','R2016b','Octave','<' ,0),...
        'HTTPS_support'    ,ifversion('>' ,0       ,'Octave','<' ,0),...
        'json'             ,ifversion('>=','R2016b','Octave','>=',7),...
        'strtrim'          ,ifversion('>=',7       ,'Octave','>=',0),...
        'accumarray'       ,ifversion('>=',7       ,'Octave','>=',0));
    checkpoint('hasFeature','CharIsUTF8')
    FeatureList.CharIsUTF8 = CharIsUTF8;
end
tf = FeatureList.(feature);
end
function v000=ifversion(v001,v002,v003,v004,v005),if nargin<2||nargout>1,...
error('incorrect number of input/output arguments'),end,persistent v006 v007 v008,if ...
isempty(v006),v008=exist('OCTAVE_VERSION','builtin');v006=[100,1] * sscanf(version,'%d.%d',2);
v007={'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;'R14SP3' 701;
'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;'R2008a' 706;'R2008b' 707;'R2009a' 708;
'R2009b' 709;'R2010a' 710;'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;'R2015b' 806;'R2016a' 900;
'R2016b' 901;'R2017a' 902;'R2017b' 903;'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;
'R2020a' 908;'R2020b' 909;'R2021a' 910;'R2021b' 911;'R2022a' 912;'R2022b' 913;'R2023a' 914};end,...
if v008,if nargin==2,warning('HJW:ifversion:NoOctaveTest',...
['No version test for Octave was provided.',char(10),...
'This function might return an unexpected outcome.']),if isnumeric(v002),v009=...
0.1*v002+0.9*ifversion_f00(v002);v009=round(100*v009);else,v010=ismember(v007(:,1),v002);if ...
sum(v010)~=1,warning('HJW:ifversion:NotInDict',...
'The requested version is not in the hard-coded list.'),v000=NaN;return,else,v009=v007{v010,2};
end,end,elseif nargin==4,[v001,v009]=deal(v003,v004);v009=0.1*v009+0.9*ifversion_f00(v009);v009=...
round(100*v009);else,[v001,v009]=deal(v004,v005);v009=0.1*v009+0.9*ifversion_f00(v009);v009=...
round(100*v009);end,else,if isnumeric(v002),v009=ifversion_f00(v002*100);if mod(v009,10)==0,...
v009=ifversion_f00(v002)*100+mod(v002,1)*10;end,else,v010=ismember(v007(:,1),v002);if ...
sum(v010)~=1,warning('HJW:ifversion:NotInDict',...
'The requested version is not in the hard-coded list.'),v000=NaN;return,else,v009=v007{v010,2};
end,end,end,switch v001,case'==',v000=v006==v009;case'<',v000=v006 < v009;case'<=',v000=v006 <=...
v009;case'>',v000=v006 > v009;case'>=',v000=v006 >=v009;end,end
function v000=ifversion_f00(v000),v000=fix(v000+eps*1e3);end
function tf=ifversion___skip_test
% Some runtimes are very twitchy about tests involving graphics. This function lists them so there
% is only a single place I need to turn them off.
persistent tf_
if isempty(tf_)
    checkpoint('ifversion___skip_test','ifversion')
    OldLinuxMatlab = isunix && ~ismac && ifversion('<','R2013a','Octave','<',0);
    checkpoint('ifversion___skip_test','ifversion')
    MacOctave = ifversion('<',0,'Octave','>',0) && ismac;
    skip = OldLinuxMatlab||MacOctave;
    
    % If the release can not be hardcoded, check two other ways to figure out whether graphics are
    % truly supported.
    if ~skip
        % If figures don't work without warnings, no graphics are likely to work.
        [str,works] = evalc(func2str(@test_figure_available)); %#ok<ASGLU>
        skip = ~works;
        if works
            % The online run tool on Matlab Answers allows figures, but doesn't allow waitbars.
            [str,works] = evalc(func2str(@test_waitbar_available)); %#ok<ASGLU>
            skip = ~works;
        end
    end
    tf_ = skip;
end
tf = tf_;
end
function tf=test_figure_available
try
    [w_msg,w_id] = lastwarn('BLANK','BLANK:BLANK');
    delete(figure);
    [w_msg,w_id] = lastwarn(w_msg,w_id);
    if strcmp(w_id,'BLANK:BLANK') && strcmp(w_msg,'BLANK')
        % No warning occurred.
        tf = true;
    else
        clc % Clear the warnings that were generated.
        error('trigger')
    end
catch
    lastwarn(w_msg,w_id); % Reset lastwarn state.
    tf = false;
end
end
function tf=test_waitbar_available
try
    delete(waitbar(0,'test if GUI is available'));
    tf = true;
catch
    tf = false;
end
end
function varargout=makedir(d)
% Wrapper function to account for old Matlab releases, where mkdir fails if the parent folder does
% not exist. This function will use the legacy syntax for those releases.
persistent IsLegacy
if isempty(IsLegacy)
    % The behavior changed after R14SP3 and before R2007b, but since the legacy syntax will still
    % work in later releases there isn't really a reason to pinpoint the exact release.
    checkpoint('makedir','ifversion')
    IsLegacy = ifversion('<','R2007b','Octave','<',0);
end
varargout = cell(1,nargout);
if IsLegacy
    [d_parent,d_target] = fileparts(d);
    [varargout{:}] = mkdir(d_parent,d_target);
else
    [varargout{:}] = mkdir(d);
end
end
function [mex_filename,fun_name]=mexname(fun_name)
%Encode runtime version information in the function name.
% This can be useful if multiple versions of Matlab or Octave need to use the
% same folder to store compiled functions, while not being compatible.
%
% This function replaces a syntax like mex_filename=[fun_name '.' mexext].
%
% Syntax:
%   mex_filename=mexname(fun_name);
%   [mex_filename,updated_fun_name]=mexname(fun_name);
persistent append
if isempty(append)
    v = version;ind=[strfind(v,'.') numel(v)];
    v = sprintf('%02d.%02d',str2double({...
        v(1:(ind(1)-1)           ) ,...
        v(  (ind(1)+1):(ind(2)-1)) }));
    v = ['v' strrep(v,'.','_')];
    if ~exist('OCTAVE_VERSION', 'builtin')
        runtime = 'MATLAB';
        type = computer;
    else
        runtime = 'OCTAVE';
        arch = computer;arch = arch(1:(min(strfind(arch,'-'))-1));
        if ispc
            if strcmp(arch,'x86_64')  ,type =  'win_64';
            elseif strcmp(arch,'i686'),type =  'win_i686';
            elseif strcmp(arch,'x86') ,type =  'win_x86';
            else                      ,type = ['win_' arch];
            end
        elseif isunix && ~ismac % Essentially this is islinux
            if strcmp(arch,'i686')      ,type =  'lnx_i686';
            elseif strcmp(arch,'x86_64'),type =  'lnx_64';
            else                        ,type = ['lnx_' arch];
            end
        elseif ismac
            if strcmp(arch,'x86_64'),type =  'mac_64';
            else                    ,type = ['mac_' arch];
            end
        end
    end
    type = strrep(strrep(type,'.',''),'-','');
    append = cell(2,1);
    append{1} = ['_' runtime '_' v '_' type];
    append{2} = [append{1} '.' mexext];
end

try % Test if fun_name is a valid name.
    if ~isvarname(fun_name),error('trigger catch block'),end
catch
    error('HJW:mexname:InvalidName',...
        'The provided input can''t be a function name')
end

mex_filename = [fun_name append{2}];
fun_name = [fun_name append{1}];
end
function OpenFileExplorer(p)
% Open the system file explorer on the indicated path.
if     ispc , system(['explorer.exe "' p '"']);
elseif ismac, system(['open "' p '"']);
elseif true , system(['xdg-open ''' p '''']);
end
end
function [opts,replaced]=parse_NameValue(default,varargin)
%Match the Name,Value pairs to the default option, attempting to autocomplete
%
% The autocomplete ignores incomplete names, case, underscores, and dashes, as long as a unique
% match can be found.
%
% The first output is a struct with the same fields as the first input, with field contents
% replaced according to the supplied options struct or Name,Value pairs.
% The second output is a cellstr containing the field names that have been set.
%
% If this fails to find a match, this will throw an error with the offending name as the message.
%
% If there are multiple occurrences of a Name, only the last Value will be returned. This is the
% same as Matlab internal functions like plot. GNU Octave also has this behavior.
%
% If a struct array is provided, only the first element will be used. An empty struct array will
% trigger an error.

switch numel(default)
    case 0
        error('parse_NameValue:MixedOrBadSyntax',...
            'Optional inputs must be entered as Name,Value pairs or as a scalar struct.')
    case 1
        % Do nothing.
    otherwise
        % If this is a struct array, explicitly select the first element.
        default=default(1);
end

% Create default output and return if no other inputs exist.
opts = default;replaced = {};
if nargin==1,return,end

% Unwind an input struct to Name,Value pairs.
try
    struct_input = numel(varargin)==1 && isa(varargin{1},'struct');
    NameValue_input = mod(numel(varargin),2)==0 && all(...
        cellfun('isclass',varargin(1:2:end),'char'  ) | ...
        cellfun('isclass',varargin(1:2:end),'string')   );
    if ~( struct_input || NameValue_input )
        error('trigger')
    end
    if nargin==2
        Names = fieldnames(varargin{1});
        Values = struct2cell(varargin{1});
    else
        % Wrap in cellstr to account for strings (this also deals with the fun(Name=Value) syntax).
        Names = cellstr(varargin(1:2:end));
        Values = varargin(2:2:end);
    end
    if ~iscellstr(Names),error('trigger');end %#ok<ISCLSTR>
catch
    % If this block errors, that is either because a missing Value with the Name,Value syntax, or
    % because the struct input is not a struct, or because an attempt was made to mix the two
    % styles. In future versions of this functions an effort might be made to handle such cases.
    error('parse_NameValue:MixedOrBadSyntax',...
        'Optional inputs must be entered as Name,Value pairs or as a scalar struct.')
end

% The fieldnames will be converted to char matrices in the section below. First an exact match is
% tried, then a case-sensitive (partial) match, then ignoring case, followed by ignoring any
% underscores, and lastly ignoring dashes.
default_Names = fieldnames(default);
Names_char    = cell(1,4);
Names_cell{1} = default_Names;
Names_cell{2} = lower(Names_cell{1});
Names_cell{3} = strrep(Names_cell{2},'_','');
Names_cell{4} = strrep(Names_cell{3},'-','');

% Allow spaces by replacing them with underscores.
Names = strrep(Names,' ','_');

% Attempt to match the names.
replaced = false(size(default_Names));
for n=1:numel(Names)
    name = Names{n};
    
    % Try a case-sensitive match.
    [match_idx,Names_char{1}] = parse_NameValue__find_match(Names_char{1},Names_cell{1},name);
    
    % Try a case-insensitive match.
    if numel(match_idx)~=1
        name = lower(name);
        [match_idx,Names_char{2}] = parse_NameValue__find_match(Names_char{2},Names_cell{2},name);
    end
    
    % Try a case-insensitive match ignoring underscores.
    if numel(match_idx)~=1
        name = strrep(name,'_','');
        [match_idx,Names_char{3}] = parse_NameValue__find_match(Names_char{3},Names_cell{3},name);
    end
    
    % Try a case-insensitive match ignoring underscores and dashes.
    if numel(match_idx)~=1
        name = strrep(name,'-','');
        [match_idx,Names_char{4}] = parse_NameValue__find_match(Names_char{4},Names_cell{4},name);
    end
    
    if numel(match_idx)~=1
        error('parse_NameValue:NonUniqueMatch',Names{n})
    end
    
    % Store the Value in the output struct and mark it as replaced.
    opts.(default_Names{match_idx}) = Values{n};
    replaced(match_idx)=true;
end
replaced = default_Names(replaced);
end
function [match_idx,Names_char]=parse_NameValue__find_match(Names_char,Names_cell,name)
% Try to match the input field to the fields of the struct.

% First attempt an exact match.
match_idx = find(ismember(Names_cell,name));
if numel(match_idx)==1,return,end

% Only spend time building the char array if this point is reached.
if isempty(Names_char),Names_char = parse_NameValue__name2char(Names_cell);end

% Since the exact match did not return a unique match, attempt to match the start of each array.
% Select the first part of the array. Since Names is provided by the user it might be too long.
tmp = Names_char(:,1:min(end,numel(name)));
if size(tmp,2)<numel(name)
    tmp = [tmp repmat(' ', size(tmp,1) , numel(name)-size(tmp,2) )];
end

% Find the number of non-matching characters on every row. The cumprod on the logical array is
% to make sure that only the starting match is considered.
non_matching = numel(name)-sum(cumprod(double(tmp==repmat(name,size(tmp,1),1)),2),2);
match_idx = find(non_matching==0);
end
function Names_char=parse_NameValue__name2char(Names_char)
% Convert a cellstr to a padded char matrix.
len = cellfun('prodofsize',Names_char);maxlen = max(len);
for n=find(len<maxlen).' % Pad with spaces where needed
    Names_char{n}((end+1):maxlen) = ' ';
end
Names_char = vertcat(Names_char{:});
end
function [opts,named_fields]=parse_print_to___get_default
% This returns the default struct for use with warning_ and error_. The second output contains all
% the possible field names that can be used with the parser.
persistent opts_ named_fields_
if isempty(opts_)
    [opts_,named_fields_] = parse_print_to___get_default_helper;
end
opts = opts_;
named_fields = named_fields_;
end
function [opts_,named_fields_]=parse_print_to___get_default_helper
default_params = struct(...
    'ShowTraceInMessage',false,...
    'WipeTraceForBuiltin',false);
opts_ = struct(...
    'params',default_params,...
    'fid',[],...
    'obj',[],...
    'fcn',struct('h',{},'data',{}),...
    'boolean',struct('con',[],'fid',false,'obj',false,'fcn',false,'IsValidated',false));
named_fields_params = fieldnames(default_params);
for n=1:numel(named_fields_params)
    named_fields_params{n} = ['option_' named_fields_params{n}];
end
named_fields_ = [...
    {'params'};
    named_fields_params;...
    {'con';'fid';'obj';'fcn'}];
for n=1:numel(named_fields_)
    named_fields_{n} = ['print_to_' named_fields_{n}];
end
named_fields_ = sort(named_fields_);
end
function opts=parse_print_to___named_fields_to_struct(named_struct)
% This function parses the named fields (print_to_con, print_to_fcn, etc) to the option struct
% syntax that warning_ and error_ expect. Any additional fields are ignored.
% Note that this function will not validate the contents after parsing and the validation flag will
% be set to false.
%
% Input struct:
% options.print_to_con=true;      % or false
% options.print_to_fid=fid;       % or []
% options.print_to_obj=h_obj;     % or []
% options.print_to_fcn=struct;    % or []
% options.print_to_params=struct; % or []
%
% Output struct:
% options.params
% options.fid
% options.obj
% options.fcn.h
% options.fcn.data
% options.boolean.con
% options.boolean.fid
% options.boolean.obj
% options.boolean.fcn
% options.boolean.IsValidated

persistent default print_to_option__field_names_in print_to_option__field_names_out
if isempty(print_to_option__field_names_in)
    % Generate the list of options that can be set by name.
    checkpoint('parse_print_to___named_fields_to_struct','parse_print_to___get_default')
    [default,print_to_option__field_names_in] = parse_print_to___get_default;
    pattern = 'print_to_option_';
    for n=numel(print_to_option__field_names_in):-1:1
        if ~strcmp(pattern,print_to_option__field_names_in{n}(1:min(end,numel(pattern))))
            print_to_option__field_names_in( n)=[];
        end
    end
    print_to_option__field_names_out = strrep(print_to_option__field_names_in,pattern,'');
end

opts = default;

if isfield(named_struct,'print_to_params')
    opts.params = named_struct.print_to_params;
else
    % There might be param fields set with ['print_to_option_' parameter_name].
    for n=1:numel(print_to_option__field_names_in)
        field_in = print_to_option__field_names_in{n};
        if isfield(named_struct,print_to_option__field_names_in{n})
            field_out = print_to_option__field_names_out{n};
            opts.params.(field_out) = named_struct.(field_in);
        end
    end
end

if isfield(named_struct,'print_to_fid'),opts.fid = named_struct.print_to_fid;end
if isfield(named_struct,'print_to_obj'),opts.obj = named_struct.print_to_obj;end
if isfield(named_struct,'print_to_fcn'),opts.fcn = named_struct.print_to_fcn;end
if isfield(named_struct,'print_to_con'),opts.boolean.con = named_struct.print_to_con;end
opts.boolean.IsValidated = false;
end
function [isValid,ME,opts]=parse_print_to___validate_struct(opts)
% This function will validate all interactions. If a third output is requested, any invalid targets
% will be removed from the struct so the remaining may still be used.
% Any failures will result in setting options.boolean.con to true.
%
% NB: Validation will be skipped if opts.boolean.IsValidated is set to true.

% Initialize some variables.
AllowFailed = nargout>=3;
ME=struct('identifier','','message','');
isValid = true;
if nargout>=3,AllowFailed = true;end

% Check to see whether the struct has already been verified.
checkpoint('parse_print_to___validate_struct','test_if_scalar_logical')
[passed,IsValidated] = test_if_scalar_logical(opts.boolean.IsValidated);
if passed && IsValidated
    return
end

% Parse the logical that determines if a warning will be printed to the command window.
% This is true by default, unless an fid, obj, or fcn is specified, which is ensured elsewhere. If
% the fid/obj/fcn turn out to be invalid, this will revert to true at the end of this function.
checkpoint('parse_print_to___validate_struct','test_if_scalar_logical')
[passed,opts.boolean.con] = test_if_scalar_logical(opts.boolean.con);
if ~passed && ~isempty(opts.boolean.con)
    ME.message = ['Invalid print_to_con parameter:',char(10),...
        'should be a scalar logical or empty double.']; %#ok<CHARTEN>
    ME.identifier = 'HJW:print_to:ValidationFailed';
    isValid = false;
    if ~AllowFailed,return,end
end

[ErrorFlag,opts.fid] = validate_fid(opts.fid);
if ErrorFlag
    ME.message = ['Invalid print_to_fid parameter:',char(10),...
        'should be a valid file identifier or 1.']; %#ok<CHARTEN>
    ME.identifier = 'HJW:print_to:ValidationFailed';
    isValid = false;
    if ~AllowFailed,return,end
end
opts.boolean.fid = ~isempty(opts.fid);

[ErrorFlag,opts.obj]=validate_obj(opts.obj);
if ErrorFlag
    ME.message = ['Invalid print_to_obj parameter:',char(10),...
        'should be a handle to an object with a writeable String property.']; %#ok<CHARTEN>
    ME.identifier = 'HJW:print_to:ValidationFailed';
    isValid = false;
    if ~AllowFailed,return,end
end
opts.boolean.obj = ~isempty(opts.obj);

[ErrorFlag,opts.fcn]=validate_fcn(opts.fcn);
if ErrorFlag
    ME.message = ['Invalid print_to_fcn parameter:',char(10),...
        'should be a struct with the h field containing a function handle,',char(10),...
        'anonymous function or inline function.']; %#ok<CHARTEN>
    ME.identifier = 'HJW:print_to:ValidationFailed';
    isValid = false;
    if ~AllowFailed,return,end
end
opts.boolean.fcn = ~isempty(opts.fcn);

[ErrorFlag,opts.params]=validate_params(opts.params);
if ErrorFlag
    ME.message = ['Invalid print_to____params parameter:',char(10),...
        'should be a scalar struct uniquely matching parameter names.']; %#ok<CHARTEN>
    ME.identifier = 'HJW:print_to:ValidationFailed';
    isValid = false;
    if ~AllowFailed,return,end
end

if isempty(opts.boolean.con)
    % Set default value.
    opts.boolean.con = ~any([opts.boolean.fid opts.boolean.obj opts.boolean.fcn]);
end

if ~isValid
    % If any error is found, enable the print to the command window to ensure output to the user.
    opts.boolean.con = true;
end

% While not all parameters may be present from the input struct, the resulting struct is as much
% validated as is possible to test automatically.
opts.boolean.IsValidated = true;
end
function [ErrorFlag,item]=validate_fid(item)
% Parse the fid. We can use ftell to determine if fprintf is going to fail.
ErrorFlag = false;
for n=numel(item):-1:1
    try position = ftell(item(n));catch,position = -1;end
    if item(n)~=1 && position==-1
        ErrorFlag = true;
        item(n)=[];
    end
end
end
function [ErrorFlag,item]=validate_obj(item)
% Parse the object handle. Retrieving from multiple objects at once works, but writing that output
% back to multiple objects doesn't work if Strings are dissimilar.
ErrorFlag = false;
for n=numel(item):-1:1
    try
        txt = get(item(n),'String'    ); % See if this triggers an error.
        set(      item(n),'String','' ); % Test if property is writable.
        set(      item(n),'String',txt); % Restore original content.
    catch
        ErrorFlag = true;
        item(n)=[];
    end
end
end
function [ErrorFlag,item]=validate_fcn(item)
% Parse the function handles. There is no convenient way to test whether the function actually
% accepts the inputs.
ErrorFlag = false;
for n=numel(item):-1:1
    if ~isa(item,'struct') || ~isfield(item,'h') ||...
            ~ismember(class(item(n).h),{'function_handle','inline'}) || numel(item(n).h)~=1
        ErrorFlag = true;
        item(n)=[];
    end
end
end
function [ErrorFlag,item]=validate_params(item)
% Fill any missing options with defaults. If the input is not a struct, this will return the
% defaults. Any fields that cause errors during parsing are ignored.
ErrorFlag = false;
persistent default_params
if isempty(default_params)
    checkpoint('parse_print_to___validate_struct','parse_print_to___get_default')
    default_params = parse_print_to___get_default;
    default_params = default_params.params;
end
if isempty(item),item=struct;end
if ~isa(item,'struct'),ErrorFlag = true;item = default_params;return,end
while true
    try MExc = []; %#ok<NASGU>
        checkpoint('parse_print_to___validate_struct','parse_NameValue')
        [item,replaced] = parse_NameValue(default_params,item);
        break
    catch MExc;if isempty(MExc),MExc = lasterror;end %#ok<LERR>
        ErrorFlag = true;
        % Remove offending field as option and retry. This will terminate, as removing all
        % fields will result in replacing the struct with the default.
        item = rmfield(item,MExc.message);
    end
end
for n=1:numel(replaced)
    p = replaced{n};
    switch p
        case 'ShowTraceInMessage'
            checkpoint('parse_print_to___validate_struct','test_if_scalar_logical')
            [passed,item.(p)] = test_if_scalar_logical(item.(p));
            if ~passed
                ErrorFlag=true;
                item.(p) = default_params.(p);
            end
        case 'WipeTraceForBuiltin'
            checkpoint('parse_print_to___validate_struct','test_if_scalar_logical')
            [passed,item.(p)] = test_if_scalar_logical(item.(p));
            if ~passed
                ErrorFlag=true;
                item.(p) = default_params.(p);
            end
    end
end
end
function [success,opts,ME,ReturnFlag,replaced]=parse_varargin_robust(default,varargin)
% This function will parse the optional input arguments. If any error occurs, it will attempt to
% parse the exception redirection parameters before returning.

% Pre-assign output.
success = false;
ReturnFlag = false;
ME = struct('identifier','','message','');
replaced = cell(0);

checkpoint('parse_varargin_robust','parse_NameValue')
try ME_ = [];[opts,replaced] = parse_NameValue(default,varargin{:}); %#ok<NASGU>
catch ME_;if isempty(ME_),ME_ = lasterror;end,ME = ME_;ReturnFlag=true;end %#ok<LERR>

if ReturnFlag
    % The normal parsing failed. We should still attempt to convert the input to a struct if it
    % isn't already, so we can attempt to parse the error redirection options.
    if isa(varargin{1},'struct')
        % Copy the input struct to this variable.
        opts = varargin{1};
    else
        % Attempt conversion from Name,Value to struct.
        try
            opts = struct(varargin{:});
        catch
            % Create an empty struct to make sure the variable exists.
            opts = struct;
        end
    end
    
    % Parse any relevant settings if possible.
    if isfield(opts,'print_to')
        print_to = opts.print_to;
    else
        checkpoint('parse_varargin_robust','parse_print_to___named_fields_to_struct')
        print_to = parse_print_to___named_fields_to_struct(opts);
    end
else
    % The normal parsing worked as expected. If print_to was provided as a field, we should use
    % that one instead of the named print_to_ options.
    if ismember('print_to',replaced)
        print_to = opts.print_to;
    else
        checkpoint('parse_varargin_robust','parse_print_to___named_fields_to_struct')
        print_to = parse_print_to___named_fields_to_struct(opts);
    end
end

% Attempt to parse the error redirection options (this generates an ME struct on fail) and validate
% the chosen parameters so we avoid errors in warning_ or error_.
checkpoint('parse_varargin_robust','parse_print_to___validate_struct')
[isValid,ME__print_to,opts.print_to] = parse_print_to___validate_struct(print_to);
if ~isValid,ME = ME__print_to;ReturnFlag = true;end
end
function [id,msg,stack,trace,no_op]=parse_warning_error_redirect_inputs(varargin)
no_op = false;
if nargin==1
    %  error_(options,msg)
    %  error_(options,ME)
    if isa(varargin{1},'struct') || isa(varargin{1},'MException')
        ME = varargin{1};
        if numel(ME)~=1
            no_op = true;
            [id,msg,stack,trace] = deal('');
            return
        end
        try
            stack = ME.stack; % Use the original call stack if possible.
            checkpoint('parse_warning_error_redirect_inputs','get_trace')
            trace = get_trace(0,stack);
        catch
            checkpoint('parse_warning_error_redirect_inputs','get_trace')
            [trace,stack] = get_trace(3);
        end
        id = ME.identifier;
        msg = ME.message;
        % This line will only appear on older releases.
        pat = 'Error using ==> ';
        if strcmp(msg(1:min(end,numel(pat))),pat)
            % Look for the first newline to strip the entire first line.
            ind = min(find(ismember(double(msg),[10 13]))); %#ok<MXFND>
            if any(double(msg(ind+1))==[10 13]),ind = ind-1;end
            msg(1:ind) = '';
        end
        pat = 'Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallbac';
        % This pattern may occur when using try error(id,msg),catch,ME=lasterror;end instead of
        % catching the MException with try error(id,msg),catch ME,end.
        % This behavior is not stable enough to robustly check for it, but it only occurs with
        % lasterror, so we can use that.
        if isa(ME,'struct') && strcmp( pat , msg(1:min(end,numel(pat))) )
            % Strip the first line (which states 'error in function (line)', instead of only msg).
            msg(1:min(find(msg==10))) = ''; %#ok<MXFND>
        end
    else
        checkpoint('parse_warning_error_redirect_inputs','get_trace')
        [trace,stack] = get_trace(3);
        [id,msg] = deal('',varargin{1});
    end
else
    checkpoint('parse_warning_error_redirect_inputs','get_trace')
    [trace,stack] = get_trace(3);
    if ~isempty(strfind(varargin{1},'%')) % The id can't contain a percent symbol.
        %  error_(options,msg,A1,...,An)
        id = '';
        A1_An = varargin(2:end);
        msg = sprintf(varargin{1},A1_An{:});
    else
        %  error_(options,id,msg)
        %  error_(options,id,msg,A1,...,An)
        id = varargin{1};
        msg = varargin{2};
        if nargin>2
            A1_An = varargin(3:end);
            msg = sprintf(msg,A1_An{:});
        end
    end
end
end
function opts=parse_warning_error_redirect_options(opts)
% The input is either:
% - an empty struct
% - the long form struct (with fields names 'print_to_')
% - the short hand struct (the print_to struct with the fields 'boolean', 'fid', etc)
%
% The returned struct will be a validated short hand struct.

if ...
        isfield(opts,'boolean') && ...
        isfield(opts.boolean,'IsValidated') && ...
        opts.boolean.IsValidated
    % Do not re-check a struct that self-reports to be validated.
    return
end

try
    % First, attempt to replace the default values with the entries in the input struct.
    % If the input is the long form struct, this will fail.
    checkpoint('parse_warning_error_redirect_options','parse_NameValue','parse_print_to___get_default')
    print_to = parse_NameValue(parse_print_to___get_default,opts);
    print_to.boolean.IsValidated = false;
catch
    % Apparently the input is the long form struct, and therefore should be parsed to the short
    % form struct, after which it can be validated.
    checkpoint('parse_warning_error_redirect_options','parse_print_to___named_fields_to_struct')
    print_to = parse_print_to___named_fields_to_struct(opts);
end

% Now we can validate the struct. Here we will ignore any invalid parameters, replacing them with
% the default settings.
checkpoint('parse_warning_error_redirect_options','parse_print_to___validate_struct')
[ignore,ignore,opts] = parse_print_to___validate_struct(print_to); %#ok<ASGLU>
end
function out=PatternReplace(in,pattern,rep)
%Functionally equivalent to strrep, but extended to more data types.
% Any input is converted to a row vector.

in = reshape(in,1,[]);
out = in;
if numel(pattern)==0 || numel(pattern)>numel(in)
    % Return input unchanged (apart from the reshape), as strrep does as well.
    return
end

L = true(size(in));
L((end-numel(pattern)+2):end) = false; % Avoid partial matches
for n=1:numel(pattern)
    % For every element of the pattern, look for matches in the data. Keep track of all possible
    % locations of a match by shifting the logical vector.
    % The last n elements should be left unchanged, to avoid false positives with a wrap-around.
    L_n = in==pattern(n);
    L_n = circshift(L_n,[0 1-n]);
    L_n(1:(n-1)) = L(1:(n-1));
    L = L & L_n;
    
    % If there are no matches left (even if n<numel(pat)), the process can be aborted.
    if ~any(L),return,end
end

if numel(rep)==0
    out(L)=[];
    return
end

% For the replacement, we will create a shadow copy with a coded char array. Non-matching values
% will be coded with a space, the first character of a match will be encoded with an asterisk, and
% trailing characters will be encoded with an underscore.
% In the next step, regexprep will be used to perform the replacement, after which indexing can be
% used to compose the final array.
if numel(pattern)>1
    checkpoint('PatternReplace','bsxfun_plus')
    idx = bsxfun_plus(find(L),reshape(1:(numel(pattern)-1),[],1));
else
    idx = find(L);
end
idx = reshape(idx,1,[]);
str = repmat(' ',1,numel(in));
str(idx) = '_';
str( L ) = '*';
NonMatchL = str==' ';

% The regular expression will take care of the lazy pattern matching. This also shifts the number
% of underscores to the length of the replacement array.
str = regexprep(str,'\*_*',['*' repmat('_',1,numel(rep)-1)]);

% We can paste in the non-matching positions. Positions where the replacement should be inserted
% may or may not be correct.
out(str==' ') = in(NonMatchL);

% Now we can paste in all the replacements.
x = strfind(str,'*');
checkpoint('PatternReplace','bsxfun_plus')
idx = bsxfun_plus(x,reshape(0:(numel(rep)-1),[],1));
idx = reshape(idx,1,[]);
out(idx) = repmat(rep,1,numel(x));

% Remove the elements beyond the range of what the resultant array should be.
out((numel(str)+1):end) = [];
end
function SelfTestFailMessage=SelfTest__error_
% Run a self-test to ensure the function works as intended.
% This is intended to test internal function that do not have stand-alone testers, or are included
% in many different functions as subfunction, which would make bug regression a larger issue.

checkpoint('SelfTest__error_','error_')
ParentFunction = 'error_';
% This flag will be reset if an error occurs, but otherwise should ensure this test function
% immediately exits in order to minimize the impact on runtime.
if nargout==1,SelfTestFailMessage='';end
persistent SelfTestFlag,if ~isempty(SelfTestFlag),return,end
SelfTestFlag = true; % Prevent infinite recursion.

test_number = 0;ErrorFlag = false;
while true,test_number=test_number+1;
    switch test_number
        case 0 % (test template)
            try ME=[];
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 1
            % Test the syntax: error_(options,msg)
            try ME=[];
                filename = tempname;
                msg = 'some error message';
                options = struct('fid',fopen(filename,'w'));
                error_(options,msg)
            catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
                fclose(options.fid);
                str = SelfTest__error_extract_message(filename);
                if ~strcmp(ME.message,msg) || ...
                        ~strcmp(str,['Error: ' msg])
                    ErrorFlag = true;break
                end
            end
        case 2
            % Test the syntax: error_(options,msg,A1,...,An)
            try ME=[];
                filename = tempname;
                msg = 'important values:\nA1=''%s''\nAn=%d';
                A1 = 'char array';An = 20;
                options = struct('fid',fopen(filename,'w'));
                error_(options,msg,A1,An)
            catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
                fclose(options.fid);
                str = SelfTest__error_extract_message(filename);
                if ~strcmp(ME.message,sprintf(msg,A1,An)) || ...
                        ~strcmp(str,sprintf(['Error: ' msg],A1,An))
                    ErrorFlag = true;break
                end
            end
        case 3
            % Test the syntax: error_(options,id,msg)
            try ME=[];
                filename = tempname;
                id = 'SelfTest:ErrorID';
                msg = 'some error message';
                options = struct('fid',fopen(filename,'w'));
                error_(options,id,msg)
            catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
                fclose(options.fid);
                str = SelfTest__error_extract_message(filename);
                if ~strcmp(ME.identifier,id) || ~strcmp(ME.message,msg) || ...
                        ~strcmp(str,['Error: ' msg])
                    ErrorFlag = true;break
                end
            end
        case 4
            % Test the syntax: error_(options,id,msg,A1,...,An)
            try ME=[];
                filename = tempname;
                id = 'SelfTest:ErrorID';
                msg = 'important values:\nA1=''%s''\nAn=%d';
                A1 = 'char array';An = 20;
                options = struct('fid',fopen(filename,'w'));
                error_(options,id,msg,A1,An)
            catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
                fclose(options.fid);
                str = SelfTest__error_extract_message(filename);
                if ~strcmp(ME.identifier,id) || ~strcmp(ME.message,sprintf(msg,A1,An)) || ...
                        ~strcmp(str,sprintf(['Error: ' msg],A1,An))
                    ErrorFlag = true;break
                end
            end
        case 5
            % Test the syntax: error_(options,ME)
            try ME=[];
                filename = tempname;
                id = 'SelfTest:ErrorID';
                msg = 'some error message';
                options = struct('fid',fopen(filename,'w'));
                try M=[];error(id,msg),catch M;if isempty(M),M=lasterror;end,end %#ok<NASGU,LERR>
                error_(options,M)
            catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
                fclose(options.fid);
                str = SelfTest__error_extract_message(filename);
                if ~strcmp(ME.identifier,id) || ~strcmp(ME.message,msg) || ...
                        ~strcmp(str,['Error: ' msg])
                    ErrorFlag = true;break
                end
            end
        case 6
            % Test the write to object option.
            % Only perform graphics-based tests on runtimes where we expect them to work.
            checkpoint('SelfTest__error_','ifversion___skip_test')
            if ifversion___skip_test,continue,end
            try ME = [];
                S.h_fig = figure('Visible','off');drawnow;
                S.h_obj = text(1,1,'test','Parent',axes('Parent',S.h_fig));
                error_(struct('obj',S.h_obj),...
                    struct(...
                    'identifier','SomeFunction:ThisIsAnIdentifier',...
                    'message',['multiline' char([13 10]) 'message']));
                close(S.h_fig)
                ErrorFlag = true;break
            catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
                close(S.h_fig)
                if ~strcmp(ME.identifier,'SomeFunction:ThisIsAnIdentifier')
                    ErrorFlag = true;break
                end
            end
        case 7
            % Test the print to function option.
            try ME = [];
                filename = [tempname '.txt'];
                fid = fopen(filename,'w');
                s_fcn = struct('h',@SelfTest__error_function_call_wrapper,...
                    'data',{{fid,'Very important error message.'}});
                error_(struct('fcn',s_fcn),...
                    struct(...
                    'identifier','SomeFunction:ThisIsAnIdentifier',...
                    'message',['multiline' char([13 10]) 'message']));
                fclose(fid);
                ErrorFlag = true;break
            catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
                fclose(fid);
                if ~strcmp(ME.identifier,'SomeFunction:ThisIsAnIdentifier')
                    ErrorFlag = true;break
                end
            end
            % Now we can test whether the contents of the file are correct.
            try
                str=SelfTest__error_extract_message(filename);
                str(str<32) = '';
                if ~strcmp(str,['Error: This <error> was caught:multilinemessageThis message ',...
                        'was included: Very important error message.'])
                    ErrorFlag = true;break
                end
            catch
                ErrorFlag = true;break
            end
        otherwise % No more tests.
            break
    end
end
if ErrorFlag
    SelfTestFlag = [];
    if isempty(ME)
        if nargout==1
            SelfTestFailMessage=sprintf('Self-validator %s failed on test %d.\n',...
                ParentFunction,test_number);
        else
            error('self-test %d failed',test_number)
        end
    else
        if nargout==1
            SelfTestFailMessage=sprintf(...
                'Self-validator %s failed on test %d.\n   ID: %s\n   msg: %s\n',...
                ParentFunction,test_number,ME.identifier,ME.message);
        else
            error('self-test %d failed\n   ID: %s\n   msg: %s',...
                test_number,ME.identifier,ME.message)
        end
    end
end
end
function SelfTest__error_function_call_wrapper(error_or_warning,ME,data)
fid = data{1};
msg = data{2};
error_(struct('fid',fid),'This <%s> was caught:\n%s\nThis message was included: %s\n',...
    error_or_warning,ME.message,msg);
end
function str=SelfTest__error_extract_message(filename)
% Extract the error message from the log file.
try
    str = fileread(filename);
catch
    str = '';return
end
ind1 = min(strfind(str,']')+2); % Strip the timestamp
ind2 = max(strfind(str,'> In')-1); % Remove the function stack.
while ismember(double(str(ind2)),[10 13 32]),ind2=ind2-1;end
str = str(ind1:ind2);
end
function SelfTestFailMessage=SelfTest__PatternReplace
% Run a self-test to ensure the function works as intended.
% This is intended to test internal function that do not have stand-alone testers, or are included
% in many different functions as subfunction, which would make bug regression a larger issue.

checkpoint('SelfTest__PatternReplace','PatternReplace')
ParentFunction = 'PatternReplace';
% This flag will be reset if an error occurs, but otherwise should ensure this test function
% immediately exits in order to minimize the impact on runtime.
if nargout==1,SelfTestFailMessage='';end
persistent SelfTestFlag,if ~isempty(SelfTestFlag),return,end
SelfTestFlag = true; % Prevent infinite recursion.

test_number = 0;ErrorFlag = false;
while true,test_number=test_number+1;
    switch test_number
        case 0 % (test template)
            try ME=[];
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 1
            try ME=[];
                x = {'abababa','aba','1'};
                expect = strrep(x{:});
                result = PatternReplace(x{:});
                if ~strcmp(expect,result),ErrorFlag = true;break,end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 2
            try ME=[];
                x = {'abababa','aba','123'};
                expect = strrep(x{:});
                result = PatternReplace(x{:});
                if ~strcmp(expect,result),ErrorFlag = true;break,end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 3
            try ME=[];
                expect = [1 4 5 3];
                result = PatternReplace([1 2 3],2,[4 5]);
                if ~isequal(expect,result),ErrorFlag = true;break,end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 4
            try ME=[];
                expect = int32([1 -10 3]);
                result = PatternReplace(int32([1 13 10 3]),int32([13 10]),int32(-10));
                if ~isequal(expect,result),ErrorFlag = true;break,end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        otherwise % No more tests.
            break
    end
end
if ErrorFlag
    SelfTestFlag = [];
    if isempty(ME)
        if nargout==1
            SelfTestFailMessage=sprintf('Self-validator %s failed on test %d.\n',...
                ParentFunction,test_number);
        else
            error('self-test %d failed',test_number)
        end
    else
        if nargout==1
            SelfTestFailMessage=sprintf(...
                'Self-validator %s failed on test %d.\n   ID: %s\n   msg: %s\n',...
                ParentFunction,test_number,ME.identifier,ME.message);
        else
            error('self-test %d failed\n   ID: %s\n   msg: %s',...
                test_number,ME.identifier,ME.message)
        end
    end
end
end
function [isLogical,val]=test_if_scalar_logical(val)
%Test if the input is a scalar logical or convertible to it.
% The char and string test are not case sensitive.
% (use the first output to trigger an input error, use the second as the parsed input)
%
%  Allowed values:
% - true or false
% - 1 or 0
% - 'on' or 'off'
% - "on" or "off"
% - matlab.lang.OnOffSwitchState.on or matlab.lang.OnOffSwitchState.off
% - 'enable' or 'disable'
% - 'enabled' or 'disabled'
persistent states
if isempty(states)
    states = {...
        true,false;...
        1,0;...
        'true','false';...
        '1','0';...
        'on','off';...
        'enable','disable';...
        'enabled','disabled'};
    % We don't need string here, as that will be converted to char.
end

% Treat this special case.
if isa(val,'matlab.lang.OnOffSwitchState')
    isLogical = true;val = logical(val);return
end

% Convert a scalar string to char and return an error state for non-scalar strings.
if isa(val,'string')
    if numel(val)~=1,isLogical = false;return
    else            ,val = char(val);
    end
end

% Convert char/string to lower case.
if isa(val,'char'),val = lower(val);end

% Loop through all possible options.
for n=1:size(states,1)
    for m=1:2
        if isequal(val,states{n,m})
            isLogical = true;
            val = states{1,m}; % This selects either true or false.
            return
        end
    end
end

% Apparently there wasn't any match, so return the error state.
isLogical = false;
end
function tf=TestFolderWritePermission(f)
%Returns true if the folder exists and allows Matlab to write files.
% An empty input will generally test the pwd.
%
% examples:
%   fn='foo.txt';if ~TestFolderWritePermission(fileparts(fn)),error('can''t write!'),end

if ~( isempty(f) || exist(f,'dir') )
    tf = false;return
end

fn = '';
while isempty(fn) || exist(fn,'file')
    % Generate a random file name, making sure not to overwrite any existing file.
    % This will try to create a file without an extension.
    checkpoint('TestFolderWritePermission','tmpname')
    [ignore,fn] = fileparts(tmpname('write_permission_test_','.txt')); %#ok<ASGLU>
    fn = fullfile(f,fn);
end
try
    % Test write permission.
    fid = fopen(fn,'w');fprintf(fid,'test');fclose(fid);
    delete(fn);
    tf = true;
catch
    % Attempt to clean up.
    if exist(fn,'file'),try delete(fn);catch,end,end
    tf = false;
end
end
function str=tmpname(StartFilenameWith,ext)
% Inject a string in the file name part returned by the tempname function.
if nargin<1,StartFilenameWith = '';end
if ~isempty(StartFilenameWith),StartFilenameWith = [StartFilenameWith '_'];end
if nargin<2,ext='';else,if ~strcmp(ext(1),'.'),ext = ['.' ext];end,end
str = tempname;
[p,f] = fileparts(str);
str = fullfile(p,[StartFilenameWith f ext]);
end
function str=unicode_to_char(unicode,encode_as_UTF16)
%Encode Unicode code points with UTF-16 on Matlab and UTF-8 on Octave.
%
% Input is either implicitly or explicitly converted to a row-vector.

checkpoint('unicode_to_char','ifversion')
persistent isOctave,if isempty(isOctave),isOctave = ifversion('<',0,'Octave','>',0);end
if nargin==1
    checkpoint('unicode_to_char','CharIsUTF8')
    encode_as_UTF16 = ~CharIsUTF8;
end
if encode_as_UTF16
    if all(unicode<65536)
        str = uint16(unicode);
        str = reshape(str,1,numel(str));%Convert explicitly to a row-vector.
    else
        % Encode as UTF-16.
        [char_list,ignore,positions] = unique(unicode); %#ok<ASGLU>
        str = cell(1,numel(unicode));
        for n=1:numel(char_list)
            checkpoint('unicode_to_char','unicode_to_UTF16')
            str_element = unicode_to_UTF16(char_list(n));
            str_element = uint16(str_element);
            str(positions==n) = {str_element};
        end
        str = cell2mat(str);
    end
    if ~isOctave
        str = char(str); % Conversion to char could trigger a conversion range error in Octave.
    end
else
    if all(unicode<128)
        str = char(unicode);
        str = reshape(str,1,numel(str));% Convert explicitly to a row-vector.
    else
        % Encode as UTF-8.
        [char_list,ignore,positions] = unique(unicode); %#ok<ASGLU>
        str = cell(1,numel(unicode)); % Create a row-vector for the result.
        for n=1:numel(char_list)
            checkpoint('unicode_to_char','unicode_to_UTF8')
            str_element = unicode_to_UTF8(char_list(n));
            str_element = uint8(str_element);
            str(positions==n) = {str_element};
        end
        str = cell2mat(str);
        str = char(str);
    end
end
end
function str=unicode_to_UTF16(unicode)
% Convert a single character to UTF-16 bytes.
%
% The value of the input is converted to binary and padded with 0 bits at the front of the string
% to fill all 'x' positions in the scheme.
% See https://en.wikipedia.org/wiki/UTF-16
%
% 1 word (U+0000 to U+D7FF and U+E000 to U+FFFF):
%  xxxxxxxx_xxxxxxxx
% 2 words (U+10000 to U+10FFFF):
%  110110xx_xxxxxxxx 110111xx_xxxxxxxx
if unicode<65536
    str = unicode;return
end
U = double(unicode)-65536; % Cast to double to avoid an error in old versions of Matlab.
U = dec2bin(U,20);
str = bin2dec(['110110' U(1:10);'110111' U(11:20)]).';
end
function str=unicode_to_UTF8(unicode)
% Convert a single character to UTF-8 bytes.
%
% The value of the input is converted to binary and padded with 0 bits at the front of the string
% to fill all 'x' positions in the scheme.
% See https://en.wikipedia.org/wiki/UTF-8
if numel(unicode)>1,error('this should only be used for single characters'),end
if unicode<128
    str = unicode;return
end
persistent pers
if isempty(pers)
    pers = struct;
    pers.limits.lower = hex2dec({'0000','0080','0800', '10000'});
    pers.limits.upper = hex2dec({'007F','07FF','FFFF','10FFFF'});
    pers.scheme{2} = '110xxxxx10xxxxxx';
    pers.scheme{2} = reshape(pers.scheme{2}.',8,2);
    pers.scheme{3} = '1110xxxx10xxxxxx10xxxxxx';
    pers.scheme{3} = reshape(pers.scheme{3}.',8,3);
    pers.scheme{4} = '11110xxx10xxxxxx10xxxxxx10xxxxxx';
    pers.scheme{4} = reshape(pers.scheme{4}.',8,4);
    for b=2:4
        pers.scheme_pos{b} = find(pers.scheme{b}=='x');
        pers.bits(b) = numel(pers.scheme_pos{b});
    end
end
bytes = find(pers.limits.lower<=unicode & unicode<=pers.limits.upper);
str = pers.scheme{bytes};
scheme_pos = pers.scheme_pos{bytes};
% Cast to double to avoid an error in old versions of Matlab.
b = dec2bin(double(unicode),pers.bits(bytes));
str(scheme_pos) = b;
str = bin2dec(str.').';
end
function [unicode,IsUTF16]=UTF16_to_unicode(UTF16)
%Convert UTF-16 to the code points stored as uint32
%
%See https://en.wikipedia.org/wiki/UTF-16
%
% 1 word (U+0000 to U+D7FF and U+E000 to U+FFFF):
%  xxxxxxxx_xxxxxxxx
% 2 words (U+10000 to U+10FFFF):
%  110110xx_xxxxxxxx 110111xx_xxxxxxxx
%
% If a second output argument is specified, an attempt is made to convert to Unicode, leaving any
% invalid encodings in place.
if nargout>=2,IsUTF16 = true;end

persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
UTF16 = uint32(UTF16);

multiword = UTF16>55295 & UTF16<57344; % 0xD7FF and 0xE000
if ~any(multiword)
    unicode = UTF16;return
end

word1 = find( UTF16>=55296 & UTF16<=56319 );
word2 = find( UTF16>=56320 & UTF16<=57343 );
try
    d = word2-word1;
    if any(d~=1) || isempty(d)
        error('trigger error')
    end
catch
    if nargout>=2
        IsUTF16 = false;
        % Remove elements that will cause problems.
        word2 = intersect(word2,word1+1);
        if isempty(word2),return,end
        word1 = word2-1;
    else
        error('input is not valid UTF-16 encoded')
    end
end

% Binary header:
%  110110xx_xxxxxxxx 110111xx_xxxxxxxx
%  00000000 01111111 11122222 22222333
%  12345678 90123456 78901234 56789012
header_bits = '110110110111';header_locs=[1:6 17:22];
multiword = UTF16([word1.' word2.']);
multiword = unique(multiword,'rows');
S2 = mat2cell(multiword,ones(size(multiword,1),1),2);
unicode = UTF16;
for n=1:numel(S2)
    bin = dec2bin(double(S2{n}))';
    
    if ~strcmp(header_bits,bin(header_locs))
        if nargout>=2
            % Set flag and continue to the next pair of words.
            IsUTF16 = false;continue
        else
            error('input is not valid UTF-16')
        end
    end
    bin(header_locs) = '';
    if ~isOctave
        S3 = uint32(bin2dec(bin  ));
    else
        S3 = uint32(bin2dec(bin.')); % Octave needs an extra transpose.
    end
    S3 = S3+65536;% 0x10000
    % Perform actual replacement.
    checkpoint('UTF16_to_unicode','PatternReplace')
    unicode = PatternReplace(unicode,S2{n},S3);
end
end
function WriteUTF8(fn,str,varargin)
%Write UTF-8 text file
%
% While it is fairly easy on new releases of Matlab to write UTF-8 files, older releases tend to
% have some difficulties creating such files.
% This function facilitates writing text and accepts inputs encoded as text (char/cellstr/string),
% as well as in a binary format (either UTF-32 or UTF-8).
%
% Syntax:
%   WriteUTF8(filename,txt)
%   WriteUTF8(filename,bin)
%   WriteUTF8(___,options)
%   WriteUTF8(___,Name,Value)
%
% Input/output arguments:
% filename:
%   A char array with either relative or absolute path.
% txt:
%   Either a char vector or matrix, or a cellstr. A string is converted to a cellstr internally.
%   Any trailing whitespace is preserved. If a cellstr is not a vector, it will be linearized.
%   Any trailing zeros will be ignored when writing the file.
% bin:
%   Either a uint32 vector/matrix, or a uint8 vector/matrix. A uint32 array will be interpreted as
%   Unicode codepoints and will be converted to UTF-8 on writing. A uint8 array will be interpreted
%   as UTF-8. There will be no check to confirm this is valid UTF-8.
%   Any trailing zeros will be ignored when writing the file.
% options:
%   A struct with the parameters. Any missing field will be filled from the defaults listed below.
%   Using incomplete parameter names or incorrect capitalization is allowed, as long as there is a
%   unique match.
%
% Name,Value parameters:
%   permission:
%      The permission used by fopen to generate a file identifier. See doc('fopen') for all
%      available options. The text mode specifier will have no effect. [default='w';]
%   EOL:
%      The end-of-line definition is usually char(10) on Unix systems and char([13 10]) on Windows
%      systems. Any char vector is allowed. [default=char(10);]
%   BOM:
%      A logical determining whether the byte order mark should be inserted. Note that this should
%      only be set to true the first time writing to a file. Writing of the BOM is skipped if ftell
%      returns a non-zero value. [default=false;]
%   SkipFileTest:
%      If set to true, this means that the validity of the filename is not checked. Since this
%      check includes a write test to see if the folder is writeable, this can represent a
%      substantial speedup for small text sizes. [default=false;]
%   print_to_con:
%      NB: An attempt is made to use this parameter for warnings or errors during input parsing.
%      A logical that controls whether warnings and other output will be printed to the command
%      window. Errors can't be turned off. [default=true;]
%      Specifying print_to_fid, print_to_obj, or print_to_fcn will change the default to false,
%      unless parsing of any of the other exception redirection options results in an error.
%   print_to_fid:
%      NB: An attempt is made to use this parameter for warnings or errors during input parsing.
%      The file identifier where console output will be printed. Errors and warnings will be
%      printed including the call stack. You can provide the fid for the command window (fid=1) to
%      print warnings as text. Errors will be printed to the specified file before being actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_obj:
%      NB: An attempt is made to use this parameter for warnings or errors during input parsing.
%      The handle to an object with a String property, e.g. an edit field in a GUI where console
%      output will be printed. Messages with newline characters (ignoring trailing newlines) will
%      be returned as a cell array. This includes warnings and errors, which will be printed
%      without the call stack. Errors will be written to the object before the error is actually
%      thrown. [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_fcn:
%      NB: An attempt is made to use this parameter for warnings or errors during input parsing.
%      A struct with a function handle, anonymous function or inline function in the 'h' field and
%      optionally additional data in the 'data' field. The function should accept three inputs: a
%      char array (either 'warning' or 'error'), a struct with the message, id, and stack, and the
%      optional additional data. The function(s) will be run before the error is actually thrown.
%      [default=[];]
%      If print_to_fid, print_to_obj, and print_to_fcn are all empty, this will have the effect of
%      suppressing every output except errors.
%      Array inputs are allowed.
%   print_to_params:
%      NB: An attempt is made to use this parameter for warnings or errors during input parsing.
%      This struct contains the optional parameters for the error_ and warning_ functions.
%      Each field can also be specified as ['print_to_option_' parameter_name]. This can be used to
%      avoid nested struct definitions.
%      ShowTraceInMessage:
%        [default=false] Show the function trace in the message section. Unlike the normal results
%        of rethrow/warning, this will not result in clickable links.
%      WipeTraceForBuiltin:
%        [default=false] Wipe the trace so the rethrow/warning only shows the error/warning message
%        itself. Note that the wiped trace contains the calling line of code (along with the
%        function name and line number), while the generated trace does not.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.0.0                                                         |%
%|  Date:    2023-07-01                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). You can see the full test matrix below.
% Compatibility considerations:
% - This is expected to work on all releases.

if nargin>=2
    [success,opts,ME] = WriteUTF8_parse_inputs(fn,str,varargin{:});
else
    success = false;
    opts.print_to = struct;
    ME = struct('identifier','HJW:WriteUTF8:nargin',...
        'message','Incorrect number of input arguments');
end
if ~success
    % If the parsing of print_to failed (which is tried first), the default will be used.
    checkpoint('WriteUTF8','error_')
    error_(opts.print_to,ME)
end

persistent BOM ExtraArgument
if isempty(BOM)
    BOM = reshape(uint8(hex2dec(['EF';'BB';'BF'])),1,3);
    checkpoint('WriteUTF8','ifversion')
    if ifversion('>=','R2006a','Octave','>=',6),ExtraArgument = {'US-ASCII'};
    else,ExtraArgument = {};end
end

% Convert the input to a cell array of uint8.
txt = WriteUTF8_prepare_binary_data( opts.str );
% Convert the EOL to uint8 and generate the BOM.
EOL = WriteUTF8_prepare_binary_data({opts.EOL});
EOL = EOL{1}; % Unwrap to avoid repeated indexing.

fid = fopen(fn,opts.permission,'n',ExtraArgument{:});
if opts.BOM && IsStartOfFile(fn,fid),fwrite(fid,BOM,'uint8');end
for n=1:numel(txt)
    fwrite(fid,txt{n},'uint8');
    fwrite(fid,EOL,'uint8');
end
fclose(fid);
end
function tf=IsStartOfFile(fn,fid)
% Check if the file is empty before writing the BOM.
tf = ftell(fid)==0;
persistent legacy
if isempty(legacy)
    checkpoint('WriteUTF8','ifversion')
    legacy = ispc && ifversion('<',7,'Octave','>',0);
end
if legacy
    % Apparently sometimes ftell is not enough, so test the file size as well.
    a = dir(fn);
    FileSize = a.bytes;
    tf = tf && FileSize==0;
end
end
function [success,opts,ME]=WriteUTF8_parse_inputs(fn,str,varargin)
% Parse the required and optional inputs to a struct. If anything fails, an ME struct will be
% created and the first output argument will be set to false.

% Attempt to match the inputs to the available options. This will return a struct with the same
% fields as the default option struct. If anything fails, an attempt will be made to parse the
% exception redirection options anyway.
checkpoint('WriteUTF8','parse_varargin_robust')
[success,opts,ME,ReturnFlag,replaced] = parse_varargin_robust(WriteUTF8_defaults,varargin{:});
if ReturnFlag,return,end

% Check optional inputs.
for k=1:numel(replaced)
    item = opts.(replaced{k});
    ME.identifier = ['HJW:WriteUTF8:incorrect_input_opt_' lower(replaced{k})];
    switch replaced{k}
        case 'permission'
            if isa(item,'string') && numel(item)==1,item = char(item);end
            % Sort the input and the allowed options to make the matching easy.
            AllowedOptions = {...
                'r' ,'w' ,'a' ,'+r' ,'+w' ,'+a' ,'A' ,'W' ,...
                'rt','tw','at','+rt','+tw','+at','At','Wt'};
            if ~isa(item,'char') || ~ismember(sort(item),AllowedOptions)
                ME.message = 'permission input parameter is not valid.';
                return
            end
            % On Windows, remove the text mode indicator to avoid adding char(13).
            if ispc,item = strrep(item,'t','');end
            opts.permission = item;
        case 'EOL'
            % This is not documented, but let's convert to char, just to be nice.
            if isa(item,'string') && numel(item)==1,item = char(item);end
            if ~isa(item,'char') || numel(item)~=length(item)
                ME.message = 'EOL input parameter is not a char vector.';
                return
            end
            opts.EOL = reshape(item,1,[]); % Just in case it wasn't a row vector.
        case 'BOM'
            checkpoint('WriteUTF8','test_if_scalar_logical')
            [isLogical,item] = test_if_scalar_logical(item);
            if ~isLogical
                ME.message = 'BOM input parameter is not a scalar logical.';
                return
            end
            opts.BOM = item;
        case 'SkipFileTest'
            checkpoint('WriteUTF8','test_if_scalar_logical')
            [isLogical,item] = test_if_scalar_logical(item);
            if ~isLogical
                ME.message = 'SkipFileTest input parameter is not a scalar logical.';
                return
            end
            opts.SkipFileTest = item;
    end
end

% Check required inputs.
if opts.SkipFileTest
    opts.fn = fn;
else % check first input
    item = fn;
    checkpoint('WriteUTF8','filename_is_valid')
    [isValid,item] = filename_is_valid(item);
    checkpoint('WriteUTF8','TestFolderWritePermission')
    if ~isValid || ~TestFolderWritePermission(fileparts(fn))
        ME.identifier = 'HJW:WriteUTF8:incorrect_input_filename';
        ME.message = 'Filename input argument is not char vector describing a writable file name.';
        return
    end
    opts.fn = item;
end
if 2==2 % check second input
    item = str;
    ME.message = 'Second input argument is not valid.';
    switch class(item)
        case 'char'
            ME.identifier = 'HJW:WriteUTF8:incorrect_input_txt';
            if ndims(item)~=2,return,end %#ok<ISMAT> This covers vectors as well.
            item(:,end+1) = 'x'; % Extend to make sure no line ends with whitespace.
            item = cellstr(item);
            for n=1:numel(item),item{n}(end) = '';end % Remove extension.
        case 'string'
            item = cellstr(item);
            item = reshape(item,1,[]);
        case 'cell'
            ME.identifier = 'HJW:WriteUTF8:incorrect_input_txt';
            if ~iscellstr(item) || ndims(item)~=2,return,end %#ok<ISCLSTR,ISMAT>
            item = reshape(item,1,[]);
        case 'uint32'
            ME.identifier = 'HJW:WriteUTF8:incorrect_input_bin';
            if ndims(item)~=2,return,end %#ok<ISMAT> This covers vectors as well.
        case 'uint8'
            ME.identifier = 'HJW:WriteUTF8:incorrect_input_bin';
            if ndims(item)~=2,return,end %#ok<ISMAT> This covers vectors as well.
        otherwise
            ME.identifier = 'HJW:WriteUTF8:incorrect_second_input';
            ME.message = 'Unsupported class for second input argument';
            return
    end
    opts.str = item;
end

% No return calls, so no errors.
success = true;
end
function opts = WriteUTF8_defaults
persistent opts_
if isempty(opts_)
    opts_ = struct(...
        'permission','w',...
        'SkipFileTest',false,...
        'EOL',char(10),...
        'BOM',false); %#ok<CHARTEN>
    % Insert the print_to option parameters.
    checkpoint('WriteUTF8','parse_print_to___get_default')
    [opts_.print_to,print_to_fieldnames_] = parse_print_to___get_default;
    for n=1:numel(print_to_fieldnames_)
        opts_.(print_to_fieldnames_{n})=[];
    end
end
opts = opts_;
end
function txt = WriteUTF8_prepare_binary_data(str)
% Convert the input to UTF-8.
if isa(str,'cell')
    checkpoint('WriteUTF8','CharIsUTF8')
    if CharIsUTF8
        % The char datatype is already UTF-8, so just convert the datatype to uint8.
        txt = str;
        for n=1:numel(str)
            row = uint8(str{n});
            while numel(row)>0 && row(end)==0,row(end)=[];end
            txt{n} = row;
        end
    else
        % Convert the char data to UTF-8.
        txt = str;
        for n=1:numel(str)
            row = str{n};
            while numel(row)>0 && row(end)==0,row(end)=[];end
            checkpoint('WriteUTF8','unicode_to_char','UTF16_to_unicode')
            if numel(row)>0,row = unicode_to_char(UTF16_to_unicode(row),false);end
            txt{n} = uint8(row);
        end
    end
elseif isa(str,'uint8') || isa(str,'uint32')
    % Convert to a cell array with UTF-8 and remove trailing zeros.
    ConvertFromUnicode = isa(str,'uint32');
    txt = cell(size(str,1),1); % Create column vector.
    for n=1:numel(txt)
        row = str(n,:);
        % Strip trailing zeros.
        while numel(row)>0 && row(end)==0,row(end)=[];end
        if numel(row)>0
            if ConvertFromUnicode
                % Convert UTF-32 to UTF-8.
                checkpoint('WriteUTF8','unicode_to_char')
                txt{n} = uint8(unicode_to_char(row,false));
            else
                txt{n} = row;
            end
        else
            txt{n} = uint8([]);
        end
    end
end
end

function out=checkpoint(caller,varargin)
% This function has limited functionality compared to the debugging version.
% (one of the differences being that this doesn't read/write to a file)
% Syntax:
%   checkpoint(caller,dependency)
%   checkpoint(caller,dependency_1,...,dependency_n)
%   checkpoint(caller,checkpoint_flag)
%   checkpoint('reset')
%   checkpoint('read')
%   checkpoint('write_only_to_file_on_read')
%   checkpoint('write_to_file_every_call')

persistent data
if isempty(data)||strcmp(caller,'reset')
    data = struct('total',0,'time',0,'callers',{{}});
end
if strcmp(caller,"read")
    out = data.time;return
end
if nargin==1,return,end
then = now;
for n=1:numel(varargin)
    data.total = data.total+1;
    data.callers = sort(unique([data.callers {caller}]));
    if ~isfield(data,varargin{n}),data.(varargin{n})=0;end
    data.(varargin{n}) = data.(varargin{n})+1;
end
data.time = data.time+ (now-then)*( 24*60*60*1e3 );
data.time = round(data.time);
end

