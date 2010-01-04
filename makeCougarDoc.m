function makeCougarDoc

%% Init variables
% rootDir = '/usr/home/rconan/Desktop/svn-rep/matlab/cougar/src';
[pathstr, name, ext, versn] = fileparts(mfilename('fullpath'));
rootDir = regexprep(pathstr,'src','');

%% Loading the function template
pathToDoc      = fullfile(rootDir,'doc','src');
pathToFile     = fullfile(pathToDoc,'class.template.html');
htmlClassLines = textread(pathToFile,'%s','whitespace','');
pathToFile     = fullfile(pathToDoc,'index.template.html');
htmlIndexLines = textread(pathToFile,'%s','whitespace','');

%% Browsing the src
files   = dir(fullfile(rootDir,'src','*.m'));
nFiles  = length(files);
nClass  = 0;
classList  = [];
for iFiles = 1:nFiles
    %     disp(files(iFiles).name(1:end-2))
    try
        s      = eval(['?',files(iFiles).name(1:end-2)]);
        nClass = nClass+1;
        classMetaData(nClass) = s;
        className = classMetaData(nClass).Name;
        classList = [classList,...
            '<li><a href="',className,'.html">',className,'</a></li>\n'];
        %     catch
        %         disp('Not a class!')
    end
end

%% Writing index files
htmlBuffer = htmlIndexLines;
htmlBuffer = regexprep(htmlBuffer,'#CLASSLIST#',classList);
writeHtmlFile(pathToDoc,'index.html',htmlBuffer{1})
 
%% Parsing meta-data
for kClass=1:nClass
    %     disp(classMetaData(kClass).Name)
    className = classMetaData(kClass).Name;
    htmlBuffer = htmlClassLines;
    htmlBuffer = regexprep(htmlBuffer,'#MAIN#',className);
    htmlBuffer = regexprep(htmlBuffer,'#CLASSNAME#',className);
    % Processing class
    lines = textread(fullfile(rootDir,'src',className),...
                '%s','delimiter','\n','whitespace','');
    k = 1;
    pattern = ['obj = ',className];
    header = [];
    while strncmp(lines{k},'%',1)
        if strncmp(lines{k}(2:end),'%',1) % Header
%             header = [header,'<p class="header">',lines{k}(3:end)];
        elseif ~isempty(strfind(lines{k},pattern))
            header = [header,'\n<p>',lines{k}(2:end)];     
        else
            header = [header,lines{k}(2:end)];
        end
        k = k+1;
    end
    header = [header,'</p>\n'];
    htmlBuffer = regexprep(htmlBuffer,'#DESCRIPTION#',header);
    % Processing propertie
    processProperties;
    %% Processing methods
    processMethods;
    %% Writting html files
    writeHtmlFile(pathToDoc,[className,'.html'],htmlBuffer{1})
end

%% Processing properties function
    function processProperties
        properties  = classMetaData(kClass).Properties;
        nProperties = length(properties);
        propertiesAttributes = [];
        for kProperties=1:nProperties
            %         disp(properties{kProperties}.Name)
            propertiesAttributes = [propertiesAttributes,...
                '<div class="attributes">\n\',...
                '<h3>',properties{kProperties}.Name,'</h3>\n',...
                '\t<ul>\n\t\t<li>GetAccess:',...
                properties{kProperties}.GetAccess,...
                '</li>\n\t\t<li>SetAccess:',...
                properties{kProperties}.SetAccess,'</li>\n'];
            if properties{kProperties}.Sealed
                propertiesAttributes = ...
                    [propertiesAttributes,'\t\t<li>Sealed</li>\n'];
            end
            if properties{kProperties}.Dependent
                propertiesAttributes = ...
                    [propertiesAttributes,'\t\t<li>Dependent</li>\n'];
            end
            if properties{kProperties}.Constant
                propertiesAttributes = ...
                    [propertiesAttributes,'\t\t<li>Constant</li>\n'];
            end
            if properties{kProperties}.Abstract
                propertiesAttributes = ...
                    [propertiesAttributes,'\t\t<li>Abstract</li>\n'];
            end
            if properties{kProperties}.Transient
                propertiesAttributes = ...
                    [propertiesAttributes,'\t\t<li>Transient</li>\n'];
            end
            if properties{kProperties}.Hidden
                propertiesAttributes = ...
                    [propertiesAttributes,'\t\t<li>Hidden</li>\n'];
            end
            if properties{kProperties}.GetObservable
                propertiesAttributes = ...
                    [propertiesAttributes,'\t\t<li>GetObservable</li>\n'];
            end
            if properties{kProperties}.SetObservable
                propertiesAttributes = ...
                    [propertiesAttributes,'\t\t<li>SetObservable</li>\n'];
            end
            propertiesAttributes = ...
                [propertiesAttributes,'\t</ul>\n</div>\n'];
            %         disp(propertiesAttributes)
        end
        htmlBuffer = regexprep(htmlBuffer,'<!--#PROPERTIES#-->',propertiesAttributes);
    end

%% Processing methods function
    function processMethods
        methods = classMetaData(kClass).Methods;
        nMethods = length(methods);
        methodsAttributes    = [];
        methodsDefiningClass = ...
            cellfun(@(x) x.DefiningClass.Name, methods,'UniformOutput',false);
        for kMethods=1:nMethods
            %         disp(methods{kMethods}.Name)
            if strcmp(methodsDefiningClass{kMethods},className)
                methodsAttributes = [methodsAttributes,...
                    '<div class="attributes">\n',...
                    '<h3>',methods{kMethods}.Name,'</h3>\n',...
                    '\t<ul>\n\t\t<li>Access:',methods{kMethods}.Access,'</li>\n'];
            else
                methodsAttributes = [methodsAttributes,...
                    '<div class="attributes">\n',...
                    '<h3>',methods{kMethods}.Name,'(',methodsDefiningClass{kMethods},')','</h3>\n',...
                    '\t<ul>\n\t\t<li>Access:',methods{kMethods}.Access,'</li>\n'];
            end
            if methods{kMethods}.Static
                methodsAttributes = ...
                    [methodsAttributes,'\t\t<li>Static</li>\n'];
            end
            if methods{kMethods}.Abstract
                methodsAttributes = ...
                    [methodsAttributes,'\t\t<li>Abstract</li>\n'];
            end
            if methods{kMethods}.Sealed
                methodsAttributes = ...
                    [methodsAttributes,'\t\t<li>Sealed</li>\n'];
            end
            if methods{kMethods}.Hidden
                methodsAttributes = ...
                    [methodsAttributes,'\t\t<li>Hidden</li>\n'];
            end
            methodsAttributes = ...
                [methodsAttributes,'\t</ul>\n</div>\n'];
            %         disp(methodsAttributes)
        end
        htmlBuffer = regexprep(htmlBuffer,'<!--#METHODS#-->',methodsAttributes);
    end

%% Writing html files function
    function writeHtmlFile(dirName,fileName,htmlLines)
        if ~isdir(dirName)
            mkdir(dirName);
        end
        dirName = fullfile(dirName,fileName);
                        disp([' ==> Writting ',dirName])
        fid = fopen(dirName,'w+');
        if fid<3
            error(['Error opening file: ',dirName]);
        end
        fprintf(fid,htmlLines);
        status = fclose(fid);
        if status<0
            error(['Cannot close ',dirName,' !'])
        end
    end

end
