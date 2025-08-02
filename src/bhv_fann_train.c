#include "fann.h"
#include "floatfann.h"
#include <unistd.h>

int main( int argc, char* argv[])
{
  const char * netfilename = argv[1];
  const unsigned int num_input = atoi(argv[3]);
  const unsigned int num_output = atoi(argv[4]);
  const unsigned int num_layers = atoi(argv[5]);;
  const unsigned int num_neurons_hidden = atoi(argv[6]);
  const float desired_error = (const float) 0.001;
  const unsigned int max_epochs = 1000;
  const unsigned int epochs_between_reports = 100;


  struct fann *ann = fann_create_standard(num_layers, num_input,
					  num_neurons_hidden, num_output);

  fann_set_activation_function_hidden(ann, FANN_SIGMOID);
  fann_set_activation_function_output(ann, FANN_SIGMOID);
  fann_set_training_algorithm        (ann, FANN_TRAIN_RPROP );
  //fann_set_train_stop_function       (ann, FANN_STOPFUNC_BIT);
  
  fann_train_on_data(ann, fann_read_train_from_file(argv[2]),
		     max_epochs, epochs_between_reports, desired_error);

  fann_save(ann, netfilename);

  fann_destroy(ann);

  return 0;
}
