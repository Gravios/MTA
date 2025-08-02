function MTAB_delete_3d_object(obj)
delete(obj.Button);
cellfun(@delete,obj.sticks);
cellfun(@delete,obj.markers);