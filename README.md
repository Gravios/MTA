MTA
===
Management and analysis package for combined Motion Capture and electrophysiology.

TODO:
@MTADrbo - rigid body objects
@MTACskeleton - collection of rbo(s)

Installation:

	1. Download the latest copy using git

       git clone https://github.com/Gravios/MTA.git 

Setup:

	2. Copy MTA/config/MTAstartup.m into your personal matlab
       directory.

        Warning: If you edit the original in the MTA/config folder, it'll cause problems
                  for others later ... well only if you're sharing it with others, if not please 
                  see the section "Sharing is caring". 

	3. Edit the paths in MTAstartup.m to correspond to your data
       locations and matlab installations:
	
        NOTE: the paths under the options of data server which are
        given as the input of MTAConfiguration, must be where you plan
        to keep your analysis data.(see my home directory for the
        example of how it will look: ~/../gravio/matlab)
        
        Note: LMU option has an additional option "project_name" which will create a

	4. Edit the DefaultArgs of your MTAstartup to conform to point to the tags of your local setup.

        Note: MTAstartup must be run each time you which to access
		data located in a different
		location. 
		
	5. Add the following line to your startup.m file in the specified
       order.

       %% MTA
       MTAstartup.m;

	6. Try it out: See MTA/analysis/example_code.m to get started.

	7. Complain to me if there are any bugs.

