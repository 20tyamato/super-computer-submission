NAME = sp01

CC = mpiicc
FLAGS = -O3 -xCORE-AVX512 -qopt-report

SOURCES = sp01.c
INCLUDES = 

OBJECTS	= $(SOURCES:.c=.o)
DEPENDS = $(SOURCES:.c=.d)

.PHONY: all clean fclean re

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

run: re
	@echo "$(NAME): Running $(NAME)..."
	@qsub sp01.bash
	@qstat

re: fclean all
