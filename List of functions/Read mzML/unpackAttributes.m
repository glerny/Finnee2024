function [allAttributes, tgtElements] = unpackAttributes(myStruct, allAttributes, elements2extract)
narginchk(1, 3)

if nargin == 1
    allAttributes = {};
    elements2extract = '';
    
elseif nargin == 2
    elements2extract = '';
    
end

fn = fieldnames(myStruct);
tgtElements = {};

for ii = 1:numel(fn)
    
    if isfield(myStruct.(fn{ii}), 'Attributes')
        
        att = fieldnames(myStruct.(fn{ii}).Attributes);
        for jj = 1:numel(att)
            
            allAttributes{end+1, 1} = [fn{ii} '/' att{jj}];
            allAttributes{end, 2} = myStruct.(fn{ii}).Attributes.(att{jj});
            
        end
    end
    
    if isfield(myStruct.(fn{ii}), 'subElements')
        
        for jj = 1:numel(myStruct.(fn{ii}).subElements)
            
            if isfield(myStruct.(fn{ii}).subElements{jj}, elements2extract)
                tgtElements{end+1} = myStruct.(fn{ii}).subElements{jj}.(elements2extract);
            end
            [allAttributes, newElements] = unpackAttributes(myStruct.(fn{ii}).subElements{jj}, allAttributes, elements2extract);
            if ~isempty(newElements)
               tgtElements = [tgtElements; newElements];
            
            end
        end
    end
end

