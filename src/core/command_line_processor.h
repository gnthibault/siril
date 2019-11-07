#ifndef SRC_CORE_COMMAND_LINE_PROCESSOR_H_
#define SRC_CORE_COMMAND_LINE_PROCESSOR_H_

gpointer execute_script(gpointer p);
int	processcommand(const char *line);
void init_completion_command();

#endif /* SRC_CORE_COMMAND_LINE_PROCESSOR_H_ */
