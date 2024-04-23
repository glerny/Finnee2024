function [s, offset] = getNodes(FID, startTag, stopTag)
%TODO: recognise opening or closing tags

s = struct();
myContent = '';

%% GET A FULL MZML ELEMENT, CORRECTING FOR NEWLINE IN TEXT FIELD AND ADDITIONAL BLANKS
offset = ftell(FID);
newLine =  strtrim(fgetl(FID));
if ~ischar(newLine), return; end

while newLine(end) ~= '>', newLine = [newLine, strtrim(fgetl(FID))]; end

if ~isempty(startTag)
    
    while 1
        IxT = strfind(newLine, ' ');
        if isempty(IxT); IxT = length(newLine); end
        if  strcmp(newLine(2:IxT-1), startTag)
            break
            
        end
        offset = ftell(FID);
        newLine =  fgetl(FID);
        if ~ischar(newLine), return; end
        newLine = strtrim(newLine);
        
        while newLine(end) ~= '>', newLine = [newLine, strtrim(fgetl(FID))]; end
    end
end

%% FIND TAG, ATTRIBUTES AND CONTENT
% check for content
if length(strfind(newLine, '>')) == 2
    IcS = strfind(newLine, '>'); IcS = IcS(1);
    IcD = strfind(newLine, '<'); IcD = IcD(2);
    myContent = newLine(IcS+1:IcD-1);
    newLine(IcS:IcD) = [];
    newLine = strrep(newLine, '/', ' /');
    newLine = strrep(newLine, '  ', ' ');
    
elseif length(strfind(newLine, '>')) > 2
    disp(newLine)
    error('Unrecognised MZML element')
    
end



IdBk = strfind(newLine, ' ');
IdEs = strfind(newLine, '="');
IdEe = strfind(newLine, '" ');
if isempty(IdBk)
    if (newLine(1:2) == "</")
        error("CheckmeMate")
        
    end
    IdBk = length(newLine);
    
    
end

if length(IdEe) < length(IdEs)
    if contains(newLine(end-1:end), '/>')
        IdEe(end+1) = length(newLine)-1;
        
    else
        IdEe(end+1) = length(newLine);
        
    end
end
currentTag_full = newLine(2:IdBk(1)-1);
currentTag = strrep(currentTag_full, '-', '_dash_');
currentTag = strrep(currentTag, ':', '_colon_');
currentTag = strrep(currentTag, '.', '_dot_');

s.(currentTag).offset = offset;

if ~isempty(myContent)
    s.(currentTag).Content = myContent;
    
end

for ii = 1:length(IdEs)
    IxS = find(IdBk < IdEs(ii), 1, 'last');
    if isempty(IxS)
        error("")
        
    end
    CurString = newLine(IdBk(IxS)+1:IdEe(ii));
    IdEq = strfind(CurString, '=');
    myAttribute = CurString(1:(IdEq(1)-1));
    myAttribute = strrep(myAttribute, '-', '_dash_');
    myAttribute = strrep(myAttribute, ':', '_colon_');
    myAttribute = strrep(myAttribute, '.', '_dot_');
    s.(currentTag).Attributes.(myAttribute) = CurString(IdEq(1)+2:end-1);
    
end

endTag = ['/', currentTag_full];
if contains(newLine, '/>') || contains(newLine, [endTag, '>'])
    
    Ix = strfind(newLine, '/>');
    if ~isempty(Ix)
        if Ix(end) == length(newLine)-1, return; end
    end
    
    Ix = strfind(newLine, [endTag, '>']);
    if ~isempty(Ix)
        if Ix(end) == length(newLine)-length(endTag)
            return
        end
    end
end

while 1
    offset = ftell(FID);
    
    newLine = fgetl(FID);
    if ~ischar(newLine), break; end
    newLine = strtrim(newLine);
    
    IxT = strfind(newLine, ' ');
    if isempty(IxT); IxT = length(newLine); end
    if  strcmp(newLine(2:IxT-1), endTag)
        break
        
    elseif strcmp(newLine(2:IxT-1), stopTag)
        fseek(FID, offset, "bof");
        break
        
    else
        fseek(FID, offset, "bof");
        
    end
    [subs, offset] = getNodes(FID, '', stopTag);
    
    %FN = fieldnames(subs);
    %if length(FN) >1, error("length field nemaes"); end
    if isfield(s, currentTag)
        if isfield(s.(currentTag), 'subElements')
            s.(currentTag).subElements{end+1} = subs;
            
        else
            s.(currentTag).subElements{1} = subs;
            
        end
    else
        s.(currentTag).subElements{1} = subs;
    end
    
    
end

