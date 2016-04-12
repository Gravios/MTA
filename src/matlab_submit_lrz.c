#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <regex.h>

#define CONF 1000

/********************************************************************

matlab_sumbit_lrz 

Type: utilitiy 
Author: Justin Graboski
Last Modified: 20160411

Purpose: submits a job to the LMU LRZ center clusters
Version: 0.01 (alpha)

Future: Add job dependeny options to reference submited jobs

*********************************************************************/


const char* program_name;
extern char** environ;

/* Prints usage information for this program to STREAM (typically 
   stdout or stderr), and exit the program with EXIT_CODE. Does not
   return.  */
void print_usage (FILE* stream, int exit_code)
{
  fprintf(stream, "\n usage: %s options [ matlab_script filebase ]", program_name);
  fprintf(stream,
          "\n Cluster Configuration:\n"
	  "        --config            filename  Name of config for job sumbition.\n"
	  "    -c  --clusters            serial  Name of cluster for job sumbition.\n"
	  "    -D  --Directory       matlab_dir  Working directory of job script.\n"
          "    -d  --dependency    option:jobid[:job_id]  option: after[any,ok,notok], singleton"
	  "    -J  --job-name        sirota_mat  Name of cluster for job sumbition.\n"
 	  "    -m  --memory              2000mb  Size of memory allocated for job.\n"
	  "    -n  --nodes                    1  Size of memory allocated for job.\n"
	  "    -p  --cpus-per-task            1  Size of memory allocated for job.\n"
          "    -r  --requeue                  0  Attempt to rerun job in case of failure.\n"
	  "    -t  --time              hh:mm:ss  Maximum wall time allocated for job.\n"        
          "\n\n Output Control:\n"
	  "    -w  --diary                       Create diary output file from matlab execution.\n"
	  "    -e  --email            EMAIL_LRZ  Default email address is EMAIL_LRZ sysenv.\n"
	  "    -o  --output            filename  Write output to file.\n"
	  "    -i  --mail-type  [all,start,end]  Email alert options.\n"
          "\n\n Matlab Controls:\n"
          "  (not functional)   --nojvm                       suppress java virtual machine at matlab startup"
          "\n\n Miscellaneous:\n"
	  "    -h  --help                        Display this usage information.\n"
          "    -l  --cluster-list                List valid clusters.\n\n"
	  "    -v  --verbose                     Print verbose messages.\n");
  exit(exit_code);
}



int main (int argc, char* argv[])
{
  program_name = argv[0];
  int next_option;
  regex_t regex;



  /* START path setup ---------------------------------------------------------*/
  char* project_dir =  getenv("PROJECT");
  char* home_dir    =  getenv("HOME");  
  char email  [40];
  snprintf(email,sizeof(email),"%s",getenv("EMAIL_LRZ"));  

  if (project_dir==NULL) { project_dir = home_dir; }

  char  sbatch_dir            [500];
  char  data_dir              [500];
  char  matlab_dir            [500];
  char  lrz_dir               [500];

  char  sbatch_base           [500];
  char  sbatch_fullpath_cmd   [500]; 
  char  sbatch_fullpath_m     [500];
  char  sbatch_fullpath_diary [500];
  char  sbatch_fullpath_out   [500];

  snprintf(data_dir,  sizeof(data_dir),   "%s/data"          ,project_dir);
  snprintf(sbatch_dir,sizeof(sbatch_dir), "%s/code/sbatch"   ,project_dir);
  snprintf(matlab_dir,sizeof(matlab_dir), "%s/code/matlab"   ,home_dir);
  snprintf(lrz_dir,   sizeof(lrz_dir),    "%s/code/lrz"      ,home_dir);

  /* END path setup -----------------------------------------------------------*/



  /* Cluster config variables */
  char job_name     [11];
  char clusters     [25];
  char node_count    [5];
  char cpus_per_task [5];
  char mem_size     [25];
  char wall_time     [15];
  char mail_type     [5];
  int  requeue  = 0;
  const char* dependency = NULL;  

  /* Misc options*/
  const char* output_filename = NULL;    // 
  int diary    = 0;                     // flag to record matlab output
  int verbose = 0;


  /* Load an optional config file */
  char conf_fullpath [500];
  snprintf(conf_fullpath,sizeof(conf_fullpath),"%s/%s",lrz_dir,"lrzc_mpp1.conf");
  FILE* fp_conf;
  char line [80];
  char flc  [80];
  char confopt [80];
  fp_conf = fopen (conf_fullpath, "r"); // open config file provided 
  while(fgets(line, 80, fp_conf) != NULL)
    {
      snprintf(flc,sizeof(flc),"%s",line);
      snprintf(confopt,sizeof(confopt),"%s",strtok(flc,"="));

      if      (!strcmp(confopt,"job_name"))      { snprintf(job_name,  sizeof(job_name), "%s",    strtok(NULL,"=")); }  
      else if (!strcmp(confopt,"clusters"))      { snprintf(clusters,  sizeof(clusters), "%s",    strtok(NULL,"=")); } 
      else if (!strcmp(confopt,"node_count"))    { snprintf(node_count,sizeof(node_count),"%s",   strtok(NULL,"=")); } 
      else if (!strcmp(confopt,"cpus_per_task")) { snprintf(cpus_per_task,sizeof(cpus_per_task),"%s",strtok(NULL,"=")); } 
      else if (!strcmp(confopt,"mem_size"))      { snprintf(mem_size,  sizeof(mem_size),  "%s",      strtok(NULL,"=")); } 
      else if (!strcmp(confopt,"wall_time"))     { snprintf(wall_time, sizeof(wall_time),"%s",    strtok(NULL,"=")); } 
      else if (!strcmp(confopt,"mail_type"))     { snprintf(mail_type, sizeof(mail_type),"%s",    strtok(NULL,"=")); } 
    }
  fclose(fp_conf);  // close the config file



  /*START of Argument Parsing */
  /* A string listing valid short options letters.  */
  const char* const short_options = "c:D:d:J:m:n:p:rt:we:o:i:hv";
  /* An array describing valid long options.  */
  const struct option long_options[] = {
    { "clusters",     1, NULL, 'c'},
    { "Directory",    1, NULL, 'D'},
    { "dependency",   1, NULL, 'd'},
    { "job-name",     1, NULL, 'J'},
    { "memory",       1, NULL, 'm'},
    { "nodes",        1, NULL, 'n'},
    { "cpus-per-task",1, NULL, 'p'},
    { "requeue",      0, NULL, 'r'},
    { "time",         1, NULL, 't'},
    { "diary",        0, NULL, 'w'},
    { "output",       1, NULL, 'o'},
    { "mail-type",    1, NULL, 'i'},
    { "help",         0, NULL, 'h'},
    { "verbose",      0, NULL, 'v'},
    { "cluster-list", 0, NULL, 'l' },
    { "config",       1, NULL, CONF},
    { NULL,      0, NULL,  0 }   /* Required at end of arry.  */
  };

  do {
    next_option = getopt_long (argc, argv, short_options,long_options, NULL);
    switch (next_option)
      {
	/* Configure cluster from file */
      case CONF: /*    --config    */      
        snprintf(conf_fullpath,sizeof(conf_fullpath),"%s/%s",lrz_dir,optarg); 
        //printf(" config: %s",conf_fullpath);
	fp_conf = fopen (conf_fullpath, "r"); 
	while(fgets(line, 80, fp_conf) != NULL)
	  {
	    snprintf(flc,sizeof(flc),"%s",line);
	    snprintf(confopt,sizeof(confopt),"%s",strtok(flc,"="));
	    if      (!strcmp(confopt,"job_name"))      { snprintf(job_name,  sizeof(job_name), "%s",    strtok(NULL,"=")); }  
	    else if (!strcmp(confopt,"clusters"))      { snprintf(clusters,  sizeof(clusters), "%s",    strtok(NULL,"=")); } 
	    else if (!strcmp(confopt,"node_count"))    { snprintf(node_count,sizeof(node_count),"%s",   strtok(NULL,"=")); } 
	    else if (!strcmp(confopt,"cpus_per_task")) { snprintf(cpus_per_task,sizeof(cpus_per_task),"%s",strtok(NULL,"=")); } 
	    else if (!strcmp(confopt,"mem_size"))      { snprintf(mem_size,  sizeof(mem_size),  "%s",      strtok(NULL,"=")); } 
	    else if (!strcmp(confopt,"wall_time"))     { snprintf(wall_time, sizeof(wall_time),"%s",    strtok(NULL,"=")); } 
	    else if (!strcmp(confopt,"mail_type"))     { snprintf(mail_type, sizeof(mail_type),"%s",    strtok(NULL,"=")); } 
	  }
	fclose(fp_conf);  // close the config file

        break;

	 /* Configure cluster from command line */
      case 'c': /* -c  --clusters  */      strcpy(clusters,     optarg);  break;
      case 'D': /* -D  --Directory */      strcpy(matlab_dir,   optarg);  break;
      case 'd': /* -d  --dependency  */    dependency = optarg;          break;
      case 'J': /* -J  --job-name  */      strcpy(job_name,     optarg);  break;
      case 'm':	/* -m  --memory    */      strcpy(mem_size,     optarg);  break;
      case 'n': /* -n  --nodes     */      strcpy(node_count,   optarg);  break;
      case 'p': /* -p  --cpus-per-task */  strcpy(cpus_per_task,optarg);  break;
      case 'r': /* -r  --requeue   */      requeue = 1;                   break;
      case 't': /* -t  --time      */      strcpy(wall_time,optarg);      break;

	/* Output Control */
      case 'w': /* -w  --diary     */      diary = 1;                     break;
      case 'e': /* -e  --email     */      strcpy(email,optarg);          break;
      case 'i': /* -i  --mail-type */      strcpy(mail_type,optarg);      break;
      case 'o': /* -o  --output    */      output_filename = optarg;      break;

	/* Miscellaneous */
      case 'h': /* -h  --help      */      print_usage (stdout, 0);
      case 'l': /*     --cluster-list */   abort();
      case 'v': /* -v  --verbose  */       verbose = 1;                   break;
      case '?': /* invalid option.  */     print_usage (stderr, 1);
      //case -1: /* Something else: unexpected.  */
        //printf ("Argument unknown code 0%o\n",next_option);
        //abort ();
      }
  }
  while (next_option != -1);

  if (verbose==1) {
    printf("\n----------------------------------------------------------------------"
           "\n Starting: %s"
           "\n----------------------------------------------------------------------"
           "\n Data path      : %s"
           "\n sbatch output  : %s"
           "\n matlab path    : %s"
           "\n lrz config path: %s"
           "\n----------------------------------------------------------------------"
	   "\n Configuration: %s\n"
 	   ,program_name,data_dir,sbatch_dir,matlab_dir,lrz_dir,conf_fullpath);
  }



  if (verbose==1) {
    printf("\n"
	   "\n job_name      : %s"
	   "\n clusters      : %s"
	   "\n node_count    : %s"
	   "\n cpus_per_task : %s"
	   "\n mem_size      : %s"
	   "\n wall_time     : %s"
	   "\n mail_type     : %s"
	   "\n requeue       : %i" 
	   "\n----------------------------------------------------------------------"
	   ,job_name,clusters,node_count,cpus_per_task,mem_size,wall_time,mail_type,requeue);
  }

  /* END of Argument Parsing */


  
  /* Create paths based on script and filebase */
  char* script_name = argv[optind];
  char  script_base     [500];
  char  script_fullpath [500];
  char* filebase = NULL;
  int isScript   = 0; 

  //printf ("optind: %i\n", optind);
  //printf ("argc: %i\n", argc);
  
  /* Designate if the first input is a script or a function
     if it is a function then use the first input as value  */
  char  aux_script_args [500];
  if (optind==argc-1) { 
    filebase = "script";
    isScript = 1; 
  }
  else {
    filebase = argv[++optind]; 
  }

  /* Concatenate further input as a list be evaluated by matlab */
  if (++optind<argc){
    do {  
      //snprintf(aux_script_args,sizeof(aux_script_args),"%s,%s",aux_script_args,argv[optind]);
      strncat(aux_script_args,",",1);
      strncat(aux_script_args,argv[optind],sizeof(argv[optind]));
      if (verbose) { printf ("\n%s\n",aux_script_args); }
    }
    while (++optind<argc);
  }
  
  char  filebase_fullpath[500];

  snprintf(script_fullpath,      sizeof(script_fullpath),      "%s/%s.m"         ,matlab_dir,script_name);
  snprintf(filebase_fullpath,    sizeof(filebase_fullpath),    "%s/%s"           ,data_dir,filebase);
  snprintf(sbatch_base,          sizeof(sbatch_base),          "%s/%s_%s"        ,sbatch_dir,script_name,filebase);
  snprintf(sbatch_fullpath_cmd,  sizeof(sbatch_fullpath_cmd),  "%s/%s_%s.cmd"    ,sbatch_base,script_name,filebase);
  snprintf(sbatch_fullpath_m,    sizeof(sbatch_fullpath_m),    "%s/%s_%s.m"      ,sbatch_base,script_name,filebase);
  snprintf(sbatch_fullpath_diary,sizeof(sbatch_fullpath_diary),"%s/%s_%s.diary"  ,sbatch_base,script_name,filebase);
  snprintf(sbatch_fullpath_out,  sizeof(sbatch_fullpath_out),  "%s/%s_%s.out"    ,sbatch_base,script_name,filebase);

  mkdir(sbatch_base,( S_IRWXU | S_IRGRP | S_IROTH ));

  if (output_filename==NULL){
    output_filename = sbatch_fullpath_out;
  }


  /* START write sbatch script ------------------------------------------------------------*/
  if (verbose==1) {
    printf("\n----------------------------------------------------------------------"
	   "\n Writing sbatch script: %s\n"
	   ,sbatch_fullpath_cmd);
  }
  
  FILE* fp_sbatch;
  fp_sbatch = fopen(sbatch_fullpath_cmd,"w");

  fprintf(fp_sbatch,"#! /bin/bash\n\n");
  fprintf(fp_sbatch,"\n");
  fprintf(fp_sbatch,
          "#SBATCH --get-user-env\n"
	  "#SBATCH -o %s\n"
	  "#SBATCH -D %s\n",
	  output_filename,matlab_dir);

  fprintf(fp_sbatch,"\n"
	  "#SBATCH -J %s\n"
	  "#SBATCH --clusters=%s\n"
	  "#SBATCH --nodes=1-%s\n"
	  "#SBATCH --cpus-per-task=%s\n"
	  "#SBATCH --mem=%s\n"
	  "#SBATCH --time=%s\n",
	  job_name,clusters,node_count,cpus_per_task,mem_size,wall_time); 


  if (dependency!=NULL){
    fprintf(fp_sbatch,"\n"
	    "#SBATCH --dependency %s\n",
            dependency);
  }

  if (!requeue) { 
    fprintf(fp_sbatch,"#SBATCH --no-requeue\n");
  } 

  fprintf(fp_sbatch,"\n"
	  "#SBATCH --mail-type=%s\n" 
	  "#SBATCH --mail-user=%s\n" 
	  "#SBATCH --export=NONE\n\n",mail_type,email);

  fprintf(fp_sbatch,
          "source /etc/profile.d/modules.sh\n"
          "module load matlab\n\n");

  fprintf(fp_sbatch,"export OMP_NUM_THREADS=%s\n\n",cpus_per_task ); 

  /* Redundant at the momemnt */
  //char pwd[500]; 
  //FILE *fp_pwd;
  //fp_pwd = popen("pwd","r");
  //fgets(pwd,sizeof(pwd),fp_pwd);
  //fclose(fp_pwd);
  //pwd[strcspn(pwd, "\r\n")] = 0;
  //fprintf(fp_sbatch,"cd %.*s/%s/\n\n",sizeof(pwd),pwd,filebase); 
  /* END Redundancy */


  fprintf(fp_sbatch,"matlab -r -nodesktop -nosplash < %s \n",sbatch_fullpath_m);

  fclose(fp_sbatch); 

  if (verbose==1) {
    printf("    ... complete"
           "\n----------------------------------------------------------------------");
  }

  /* END write sbatch script --------------------------------------------------------------*/




  /* START write matlab script ------------------------------------------------------------*/
  if (verbose==1) {
    printf("\n----------------------------------------------------------------------"
	   "\n Writing matlab script: %s\n"
	   ,sbatch_fullpath_m);
  }

  FILE* fp_mat;
  fp_mat = fopen(sbatch_fullpath_m,"w");

  fprintf(fp_mat,"\ndisp([\'Executing batch with file : ', which(\'%s\'),\'\\n\']);\n\n",script_name);
  fprintf(fp_mat,"\naddpath %s\n",matlab_dir);
  //fprintf(fp_mat,"cd %s/%s\n",data_dir,filebase);
  if (diary) { fprintf(fp_mat, "diary %s\n\n",sbatch_fullpath_diary); }


  if (isScript) {
    fprintf(fp_mat,"%s\n",script_name);
  }
  else{
    fprintf(fp_mat,"%s(\'%s\'%s)\n",script_name,filebase,aux_script_args);
  }

  if (diary) {fprintf(fp_mat, "diary off\n",sbatch_fullpath_diary); }

  fprintf(fp_mat, "quit\n",sbatch_fullpath_diary);
  fclose(fp_mat);

  if (verbose==1) {
    printf("    ... complete"
           "\n----------------------------------------------------------------------");
  }

  /* END write matlab script --------------------------------------------------------------*/




  /* START sbatch submission */
  if (verbose==1) {
    printf("\n----------------------------------------------------------------------"
	   "\n Submitting script to cluster: %s\n"
	   ,clusters);
  }

  char cmd_start[500]; 
  char cmd_out  [500];
  char *cmd_jid = cmd_out;
  FILE *fp_sysout;

  snprintf(cmd_start,sizeof(cmd_start),"sbatch %s",sbatch_fullpath_cmd); 
  fp_sysout = popen(cmd_start,"r");
  fgets(cmd_out,sizeof(cmd_out),fp_sysout);
  fclose(fp_sysout);

  if (verbose==1) {
    printf("    ... complete"
           "\n----------------------------------------------------------------------");
  }
  /* END sbatch submission */


  /* Log job in active list */
  int reti;
  size_t nmatch = 1;
  char job_id  [20];
  regmatch_t pmatch[1];
  reti = regcomp(&regex,"job \\([[:digit:]]*\\) on",0);
  reti = regexec(&regex,cmd_out,nmatch,pmatch,0);
  regfree(&regex);
  snprintf(job_id, sizeof(job_id),"%.*s",
           pmatch[0].rm_eo-pmatch[0].rm_so-7,cmd_jid+pmatch[0].rm_so+4);
  if (verbose==1) {
    printf("\n %s "
	   "\n Job Id : %s\n"
	   "\n Expect email updates at: %s \n\n"
	   ,cmd_out,job_id,email);
  }
  printf("%s", job_id);

  return 0;
}
