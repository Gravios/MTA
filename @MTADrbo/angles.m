function ang = angles(Data);
ang = quaternion2rad(Data(:,:,5:8));