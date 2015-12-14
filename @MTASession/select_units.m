
function units = select_units(Session,query) %#ok<INUSL>
%units = select_units(Session,query)
%Return unit ids based on a query   
%
%  query: cellArray(formated msg) 
%
%  query example:
% 
%    query = {{PlaceField_walk,'maxRate',@gt,3},@or,{PlaceField_rear,'maxRate',@gt,3}}

unit_status = [];
fh_list = {};
% If the query is a single criteria wrap in a cell
if ~iscell(query{1})
    query = {query};
end
% parse query
while numel(query)~=0,
    if iscell(query{1}),
        obj = query{1}{1};                
        objmeta = metaclass(obj);

        if numel(query{1})>1,
            prop = query{1}{2};
            operator = query{1}{3};
            selection_value = query{1}{4};
        end
    else
        obj = query{1};                
        objmeta = metaclass(obj);
    end

    switch objmeta.Name
      case 'MTAApfs'
        % only really meant for the property max rate
        % at the moment
        target_prop = obj.data.(prop);
        value = cellfun(@max,target_prop,'UniformOutput',false);
        value(cellfun(@isempty,value))={[0]};
        value = cell2mat(value);                    
        unit_status = cat(2,unit_status,operator(value,selection_value));
      case 'struct'
        value = obj.(prop)
        unit_status = cat(2,unit_status,operator(value(:),selection_value));
      case 'function_handle'
        % stores the logical operators in a cell array
        fh_list{end+1} = obj;
    end
    query(1) = [];
end

while size(unit_status,2)>1,
    unit_status(:,2) = fh_list{1}(unit_status(:,1),unit_status(:,2));
    unit_status(:,1) = [];
    fh_list(1) = [];
end

units = find(unit_status);
end
