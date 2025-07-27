NAME = sp05

CC = mpiicc
FLAGS = -O3 -xCORE-AVX512 -qopt-report -qopt-report-phase=vec

SOURCES = sp05.c
INCLUDES = 

OBJECTS	= $(SOURCES:.c=.o)
DEPENDS = $(SOURCES:.c=.d)

.PHONY: all clean fclean re run submit

all: $(NAME)

$(NAME): $(OBJECTS)
	@$(CC) $(FLAGS) $(OBJECTS) -o $(NAME)
	@echo "$(NAME): object files were created"
	@echo "$(NAME): $(NAME) was created"

-include $(DEPENDS)

%.o : %.c
	@$(CC) $(FLAGS) -c $< -o $@
	@echo ".\c"

clean:
	@rm -rf $(OBJECTS) $(DEPENDS)
	@echo "$(NAME): object and dependency files were deleted"

fclean: clean
	@rm -f $(NAME)
	@echo "$(NAME): $(NAME) was deleted"

submit: re
	@echo "$(NAME): Submitting job to queue..."
	@qsub sp05.bash
	@echo "$(NAME): Job submitted. Checking queue status..."
	@sleep 2
	@qstat

run: submit

re: fclean all
