function index = get_pose_index(Data)
index = round((Data.poseIndex-Data.origin)*Data.sampleRate);