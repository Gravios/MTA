function MTAB_delete_3d_object(obj)
cellfun(@delete,obj.sticks)
cellfun(@delete,obj.markers)