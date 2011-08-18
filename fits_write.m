function fits_write(filename, data, varargin)
global BITPIX transposeRC
%
% mfits: A simple fits package fo MATLAB.
%
% Author: Jason D. Fiege, University of Manitoba, 2010
% fiege@physics.umanitoba.ca
%
% The main component is a FITS writer fits_write to complement MATLAB's
% built-in fits reader.  fits_read is just a wrapper for fitsread.  fits_info
% returns the header information in a more convenient form than the
% built-in fitsinfo function.  Additional functionality is provides by
% getFitsValue and setFitsValue to get and set a single header field.
%
% -------------------------
% Usage:
% -----
% 1st argument: filename: string containing the name of the file to write.  '.fits' will
% be appended if the file name does not end in '.fits' or '.FITS'.
%
% 2nd argument: data: an N-dimensional array, usually of 32-bit floats or 64-bit doubles.
% Complex data types are not currently supported, but this may be added in
% the the near future.
%
% 3rd argument: header (optional): A fits header can be optionally supplied
% in one of 3 formats.
%
% Header format #1: The header is a structure loaded directly from MATLAB's
% fitsinfo command: header=fitsinfo('myFITS.fits');  In this case, the
% keyword/values are held in a cell array: header.PrimaryData.Keywords.  This
% cell array is extracted and header format #2 is used.  For example, if
% 'myFITS.fits' is a fits file, then the following commands read it in and
% re-write identical data with a copy of the same header:
%
% data=fitsread('myFITS.fits');
% header=fitsinfo('myFITS.fits');
% fitswrite('myFITS_copy.fits', data, header);
%
% Header format #2: The header may be a cell array, in the same format
% as the headers returned by MATLAB's fitsinfo command in the
% header.PrimaryData.Keywords field.  For example, if 'myFITS.fits' is a fits
% file, then the following commands read it in and re-write identical data
% with a copy of the same header:
%
% data=fitsread('myFITS.fits');
% header=fitsinfo('myFITS.fits');
% fitswrite('myFITS_copy.fits', data, header.PrimaryData.Keywords);
%
% This format preserves all bare FITS keywords (with no value field) and
% all comments.
%
% Header format #3: The header may be a MATLAB structure in the format
% header.keyword1=value1;
% header.keyword2=value2;
% header.keyword3=value3;
% 
% This format is compatible with the MFITSIO package.  Note that it does
% not permit bare keywords or comment fields.
%
% Certain mandatory FITS fields may be added of they are not present in
% the supplied header.  Please see the function 'requiredKeywords' in this
% file for details.  These required keywords will be used to generate a
% minimal FITS header if the optional header field is not supplied.
%
% ==============================
% Modify these lines to suit your data and machine architecture.
%
% FITS file BITPIX field.  Common choices are -32 for single precision
% data, and -64 for double precision.
BITPIX=-32; % [-32, -64]
%
% machineFormat describes the "endian-ness" of your machine.  Valid choices
% are 'native', 'L', 'B', and others described in MATLAB's documentation
% for 'fopen'.  Type 'help fopen' on the command line for details.
% *** N.B. NOTE THE MATLAB'S FITSREAD COMMAND ALWAYS USES machineFormat='B'.
% THIS OPTION IS THEREFORE REQUIRED FOR COMPATIBILITY WITH FITSREAD.M.
% USE WITH OTHER FITSREADERS MAY REQUIRE OTHER OPTIONS. ***
machineFormat='B'; % ['native', 'L', 'B']
%
% MATLAB's built-in FITS reader seems to transpose the first 2 dimensions
% (Rows and Columns) of a multi-dimensional data set.  Set 'transposeRC=true'
% to write files that are compatible with this behaviour, or 'transposeRC=false' if
% you do not wish to transpose.  If your data seems to have the X and Y
% axis transposed when viewed with an external viewer, then you have
% probably made the wrong choice.
transposeRC=true;
% ==============================

% Retrieve or initialize the header.
if isempty(varargin)
    header=[];
else
    header=varargin{1};
end

% Check if the header is in format #1.  Extract the required keyword/value
% cell array if it is.
if isstruct(header) && isfield(header, 'PrimaryData') &&...
        isfield(header.PrimaryData, 'Keywords')
    header=header.PrimaryData.Keywords;
end

% Obtain the size of the data and transpose rows and columns for
% compatibility with fitsread if requested.
szData=size(data);
if length(szData) > 1 && transposeRC
    order=1:length(szData);
    order(1:2)=[2,1];
    data=permute(data, order);
end

% Generate an N x 80 character array hStr containing the FITS header
% information.  'struct2str' is called if the supplied header is a
% structure, or 'cell2str' if it is in cell array format.
if isstruct(header)
    [hStr, hReq]=struct2str(header, szData);
elseif iscell(header)
    [hStr, hReq]=cell2str(header, szData);
else
    % No header was supplied or the format is wrong.
    hStr=[];
    hReq=requiredKeywords(szData);
end
% Mandatory header fields *that were not over-written by the supplied
% header* are returned in the hReq.  They are prepended to the FITS
% character array below.
%
% Add the mandatory 'END' line.
hStr=addLine(hStr, 'END');

% Add any remaining mandatory fields that have not already been
% over-written by the supplied header.
if ~isempty(hReq)
    hStr=[struct2str(hReq, szData); hStr];
end

% Determine the 'precision' string that corresponds to the BITPIX field.
prec=getPrecision;

% Append the file extension '.fits' if the supplied name does not already
% end in '.fits' or '.FITS'.
if isempty(regexpi(filename, '\.fits$'))
    filename=[filename, '.fits'];
end
fid=fopen(filename, 'w+', machineFormat); % Open file for writing.
fwrite(fid, hStr', 'char'); % Write the header.
%
% Pad the header with 8-bit zeros to satisfy the FITS format requirement
% that data are written to 2880 bye blocks.
padFile(fid, 8, machineFormat);
%
fwrite(fid, data, prec, 0, machineFormat); % Write the image data.
%
% Pad the image with zeros (of the same BITPIX format) to satisfy the FITS
% format requirement that data are written to 2880 bye blocks.
padFile(fid, BITPIX, machineFormat)
fclose(fid);

% ------------------------------

function padFile(fid, BITPIX, machineFormat)
% Pad the file with zeros - FITS standard requires block sizes of 2880 bytes.
modBlock=mod(ftell(fid),2880);
if modBlock > 0
    NPad=(2880-modBlock)/(abs(BITPIX)/8);
    padding=zeros(1,NPad);
    fwrite(fid, padding, getPrecision(BITPIX), 0, machineFormat);
end

% ------------------------------

function  [hStr, hReq]=cell2str(header, szData)
% Generates a FITS header from a cell array in the same format as returned
% by MATLAB's fitsinfo function.  hReq is a list of *unused* mandatory fields
% in the FITS standard.
[nr,nc]=size(header);
hReq=requiredKeywords(szData);
hStr='';
for i=1:nr
    keyword=header{i,1};
    value=header{i,2};
    %
    if nc > 2 && ~isempty(header{i,3})
        comment=strtrim(char(header{i,3}));
    else
        comment='';
    end
    %
    % Generate a single header line and append it to hStr.  Note that the
    % mandatory FITS fields in hReq are removed as they are used.  This
    % will only contain unused fields at the end.
    [hLine, hReq]=makeHeaderLine(keyword, value, comment, hReq);
    hStr=[hStr; hLine]; %#ok
end

% ------------------------------

function  [hStr, hReq]=struct2str(header, szData)
% Generates a FITS header from a MATLAB structure in the same format
% used by mfitsio.  hReq is a list of *unused* mandatory fields
% in the FITS standard.
keywords=fields(header);
nr=length(keywords);
% No comments or bare keywords are possible when using structures.
hReq=requiredKeywords(szData);
hStr='';
for i=1:nr
    keyword=keywords{i};
    value=header.(keyword);
    %
    % Generate a single header line and append it to hStr.  Note that the
    % mandatory FITS fields in hReq are removed as they are used.  This
    % will only contain unused fields at the end.
    [hLine, hReq]=makeHeaderLine(keyword, value, '', hReq);
    hStr=[hStr; hLine]; %#ok
end

% ------------------------------

function [hLine, hReq]=makeHeaderLine(keyword, value, comment, hReq)
% Generate a single line of the FITS header.  hReq is a structure containing the
% list of mandatory keywords in the FITS standard.  Fields are removed from
% hReq as they are encountered and over-written using the supplied header
% information.
%
% Get the keyword and check that is is not 'END'.  The mandatory 'END'
% keyword is added automatically.  The ensures that the 'END' keyword is
% not duplicated.
keyword=upper(strtrim(char(keyword)));
if strcmpi(keyword, 'END')
    hLine=[];
    return
end
%
% Trim the keyword to the allowed 8 bytes if it is too long.
lenKey=length(keyword);
if lenKey > 8
    keyword=keyword(1:8);
    lenKey=8;
end
%
% Check whether keyword is in the list of required keywords.  If it is,
% then remove the keyword from hReq.
if isfield(hReq, keyword)
    value=hReq.(keyword);
    hReq=rmfield(hReq, keyword);
end
%
% Validate the value field.  This converts logicals to 'T' or 'F', and
% converts numeric data to strings.  The FITS standard allows a maximum of
% 70 characters for the value field.  The value is truncated if this is
% exceeded.
value=validate(value);
lenValue=length(value);
if lenValue > 70
    value=value(1:70);
    lenValue=70;
end
%
% If there is a comment field supplied, make sure that it is prepended by a
% space and a '/' character.  The comment is truncated if the entire line
% extends past 80 characters.
comment=strtrim(comment);
if any(comment) && comment(1) ~= '/'
    comment=[' /', comment];
end
lenComment=min(length(comment), 70-lenValue);
comment=comment(1:lenComment);
%
hLine=repmat(' ', 1, 80);
hLine(1:lenKey)=keyword;
if lenValue > 0
    hLine(9:10)='= ';
    hLine(11:11+lenValue-1)=value;
end
if lenComment > 0
    hLine(11+lenValue:10+lenValue+lenComment)=comment;
end

% ------------------------------

function value=validate(value)
% Formats the 'value' field.
if isempty(value)
    value='';
elseif ischar(value)
    % Make sure that the value string is surrounded by single quotes.
    if value(1) == ''''
        value(1)=[];
    elseif value(end) == ''''
        value(end)=[];
    end
    value=['''', strtrim(value), ''''];
elseif islogical(value)
    % Convert logicals to 'T' or 'F'.
    if value
        value='T';
    else
        value='F';
    end
elseif isnumeric(value)
    % Convert numerical values to a string representation.
    value=num2str(value);
end
value=upper(strtrim(value));

% ------------------------------

function hStr=addLine(hStr, hLine)
% Just adds a line to the header character array.
%
lenLine=length(hLine);
if lenLine > 80
    hLine=hLine(1:80);
    lenLine=80;
end
endLine=repmat(' ', 1, 80);
endLine(1:lenLine)=hLine;
hStr=[hStr; endLine];

% ------------------------------

function prec=getPrecision(varargin)
global BITPIX
% Determine the precision string from the BITPIX field.  This is used by
% fopen and fwrite.
if isempty(varargin)
    BPX=BITPIX;
else
    BPX=varargin{1};
end
%
switch BPX
    case 8
        prec='uint8';
    case 16
        prec='int16';
    case 32
        prec='int32'; %  'integer*4'
    case 64
        prec='int64'; % 'integer*8'
    case -32
        prec='single'; % 'float32'; % real*4
    case -64
        prec='double'; % 'float64'; % real *8
    otherwise
        disp('Invalid BITPIX value.');
        BITPIX=-32;
        prec='single'; % 'float32'; % real*4
end

% ------------------------------

function h=requiredKeywords(szData)
global BITPIX transposeRC
% The list of mandatory keywords specified by the FITS standard.
%
h.SIMPLE=true;
h.BITPIX=BITPIX;
h.NAXIS=length(szData);
%
% Swap axes 1 & 2 if the axes are transposed.
if h.NAXIS > 1 && transposeRC
    szData([2,1])=szData([1,2]);
end
%
% Generate the keywords for NAXIS1, NAXIS2, etc.
for i=1:h.NAXIS
    h.(['NAXIS', num2str(i)])=szData(i);
end
