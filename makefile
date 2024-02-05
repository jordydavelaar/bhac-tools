CFLAGS = -openmp  -O2   -Wall -Wno-unused-but-set-variable 
LDFLAGS = -openmp -lm -lfftw3 -I/opt/homebrew/Cellar/gsl/2.7.1/include -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas  

OBJDIR=build

TARGET=BHAC-tools

SOURCES=main.c model.c metric.c GRmath.c constants.c uniform.c
OBJECTS := $(patsubst %.c,$(OBJDIR)/%.o,$(SOURCES))

all: create_directories $(SOURCES) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

create_directories:
	@test -d $(OBJDIR) || mkdir -v $(OBJDIR)


$(OBJECTS): $(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS)  -c $^ -o $@


clean:
	rm -rf $(OBJECTS) $(TARGET)
