#ifndef SRC_CORE_COMMAND_LINE_PROCESSOR_H_
#define SRC_CORE_COMMAND_LINE_PROCESSOR_H_

gpointer execute_script(gpointer p);
int	processcommand(const char *line);
void init_completion_command();
sequence *load_sequence(const char *name, char **get_filename);

#endif /* SRC_CORE_COMMAND_LINE_PROCESSOR_H_ */
