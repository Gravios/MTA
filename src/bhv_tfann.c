#include "/usr/local/include/fann.h"
#include <unistd.h>

int main( int argc, char* argv[])
{
  const unsigned int num_input = atoi(argv[1]);
  const unsigned int num_output = atoi(argv[2]);
  const unsigned int num_layers = atoi(argv[3]);;
  const unsigned int num_neurons_hidden = atoi(argv[4]);
  const float desired_error = (const float) 0.001;
  const unsigned int max_epochs = 5000;
  const unsigned int epochs_between_reports = 1000;

  struct fann *ann = fann_create_standard(num_layers, num_input,
					  num_neurons_hidden, num_output);

  fann_set_activation_function_hidden(ann, FANN_SIGMOID_SYMMETRIC);
  fann_set_activation_function_output(ann, FANN_SIGMOID_SYMMETRIC);

  fann_train_on_file(ann, "xor.data", max_epochs,
		     epochs_between_reports, desired_error);

  fann_save(ann, "xor_float.net");

  fann_destroy(ann);

  return 0;
}
