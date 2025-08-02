function [in, out, err] = popen(cmd)
rt = java.lang.Runtime.getRuntime();
p = rt.exec(cmd);
in  = java.io.BufferedReader(java.io.InputStreamReader(p.getInputStream()));
err = java.io.BufferedReader(java.io.InputStreamReader(p.getErrorStream()));
out = java.io.PrintWriter(p.getOutputStream());
