MTA
===
is a:
	1. Magic Toolbox of Analysis
	2. Toolbox created for analyzing data. 
	3. Management and analysis package for combined Motion Capture and electrophysiology.

TODO:
@MTADrbo - rigid body objects
@MTACskeleton - collection of rbo(s)

TODO MAJOR:
change zero filling to nan filling where data doesn't exist



Warning:

	Please consult your doctor before applying this code to your
    data. Side effects may include confusion, nausea, restless leg
    syndrome and occational bouts of procrastination. If you
    experience any absent mindedness, do not painic this is normal.



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
		
	5. Add the following lines to your startup.m file in the specified
       order.

       Note: The order is important. Some functions within MTA must
       override their counterpart in the labbox.

       %% Labbox
       addpath(genpath('/gpfs01/sirota/homes/share/matlab/labbox'));
       rmpath(genpath('/gpfs01/sirota/homes/share/matlab/labbox/.git'));

       %% MTA
       MTAstartup.m;

	6. Try it out: See MTA/analysis/example_code.m to get started.

	7. Complain to me if there are any bugs, but only if it has 7 or
       11 legs, if it has 6 use your entomology skills to identify it.


Sharing is Caring:

	1. Sharing is caring
	2. Those who share should not be careless
	3. If you don't care then you're probably the reason we can't have
       nice things.



The Nature of Behavior:

	1. How should we think about time scales of behavior? 
	   
	   Cyclical vs Phasic
	   Reflexes ( short, quick, atomic, response to external environment )
	   Neuroendocrine ( modulatory, contineous signal of internal state )
       Developmental ( age dependent alteration of ion channel subunit
                       expression; hormonal shifts during sexual
                       maturation)
       Social ( Interaction of multiple organisms )
	   
  
     What is a sequence of behavior versus an individual behavior? 
	   The study of behavior should probably have a structural
	   system. Adapting the structural schema from protein folding:
     	   Primary: CPG's and Phasic neural responses
		   Secondary: primary sequences group into neural state
			   dependent structures ( e.g. Sequence of rythmic
  					  head movements used to gather environmental
				      information; gather/gaurd food; copulate;
	                  self-defence)
	       Tertiary: secondary sequences group into Behavioral pattern
			        (i.e. Grooming: Forelimb lick to head groom to ear scratch)
		   Quanternery: tertiary squences group into Behavioral
			           strategies (e.g. move forward; sweeping sensory
					   apparatuses through environment; observation,
                       computation, reaction -> consume resource; move left)
	 
	 
	 How do we motivate the choice of timescales? 
	   Observation of an organisims habitat would contain the most
	   relavent information regarding timescale. 
	   
	   Many vertibrate organisms have specific brain state sequences
	   predetermined by a combination of genetic and neuroendocrine
	   clocks. All biological clocks are an evolutionary heritage of
	   enivironmental pressures. Cycles within the habitat
	   (e.g. night/day, high/low tide, seasons, preditor/prey
	   equilibrium,...) have molded the network connectivity to be
	   optimally responsive to the physical/chemical tendencies of
	   their historical environments. Each animals "natural" habitat
	   should contain clues to the optimal timescales of observation.

  	   Determine the minimal sampling schema for detecting behavioral
	   developmental transitions, and the interaction of developmental,
  	   environmental and interoceptive contexts.
	   
	 
	 Can neural recordings help inform the choice of timescales?
       Identification of the central patern generators of differnte
       behaviors/motivatiors timescales (e.g. single cycle of walking;
       circadian rhythm; menstral cycle; migration of a
       species). "Neural recording" is vauge and stepping into more
       complicated sesnsors. A list of possible sensors types could be:
           Electrophisiology ( Oscillations of network acivity and
                               single units)
		   Neuromodulatory ( Attention, stress, social programs  )
		   Neuroendocrine ( homeostatic, sexual drive, stress, ...)
           Neurogenetic ( developmental programs, genetic clocks,...)
	   Neural recordings could span seconds to days to years or even
       decades depending on the behavior or modification of the
       organisims of intrest.

      
	  
  2. Are there “atoms” of behavior?
       Yes. (central pattern generators have been around for quite some
       time now).
  
  3. Variability: What does it mean when quantifying behavior? 
       Though atomic behaviors exist their interation may not have a
       locked phase and also employment/magnitude in/of response to environmental
       change may be continuous.
 
     How do we think about variability versus noise versus stereotypy? 
       Stereotype ... do you mean a statistical mean?
	   
	 Can you have a behavior without stereotypy? 
       Can you have a scientific study without statistics (stereotypes)?
	 
	 If an organism does something only once but never again, is this a “behavior”?
       If each organism of a species has the same/stereotypical
       behavioral response to a context then yes.

       All organisms only die once. 	  	 
	   
       Yes, Many invertabrates can attest to this... well not really
	   since their mates bit their heads off. 

       Most of the western world has been circumsized, I'm pretty sure
       they only do that once :)


  4. Periodicity: Can an aperiodic sequence be a behavior?
       Yes. see response to question 3


  5. Is behavior discrete or continuous?
       ... umm, I'm not going to waste my time with this one.


  6. Can behavioral states be enumerated or is the best we can do to
     describe some scaling law that relates the time scale of observation
     to the number of states (e.g. Coastline of Britain)?
       ... umm, again, I'm not going to waste my time with this one.


 Using Behavior to Study the Brain:
 
  1. To what extent do we expect the nature of behavior to reflect the
     nature of the brain? E.g. do discrete behaviors imply discrete brain
     states? 
	 
	 Should the dimensionality of behavior match dimensionality
     of neural dynamics?
	 
	 
  2. Behavior time scales versus neural time scales; how do we think
     about differences in these scales? Should long-lived behavioral
     states necessarily have long-lived neural correlates? Similarly,
     should high-frequency behaviors (flapping wings, for example) have
     correlates at those frequencies?
	 
	 
  3. What are the paths forward for relating neural activity to unbiased
     measures of behavior, and what are the obstacles?


 Placing unbiased studies of Behavior in Context:
  
  It is all subjective.
 
  1. What has unbiased quantification of behavior taught us that other
     methods could not have?
	   The qualities of behevaior we choose to quantify may or more
       often do not contain the maximum information regarding the
       mechanisms of generation. The selection and driving of CPG's
       and phasic neural responses to simple stimuli has been very
       successfully analysed in simplified contexts, but a large
       amount of variability in behaviors.
	 
  2. What do we hope it will teach us in the future about both behavior
     and the brain
         


