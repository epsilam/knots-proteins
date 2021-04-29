CXX      := g++
CXXFLAGS := -Wall -std=c++2a
LIBS     := -pthread
SRCDIR   := src
OBJDIR   := obj
BINDIR   := bin
TARGET   := $(BINDIR)/main

# Filename extension of source files
SRCEXT   := cc

SOURCES  := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS  := $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SOURCES))

# Link
$(TARGET): $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(CXX) $^ -o $(TARGET) $(LIBS)

# Compile source files
$(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) 
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -rf $(BINDIR) $(OBJDIR) $(ERRLOG)

.PHONY: cleanvid
cleanvid:
	rm -rf $(VIDDIR)

run: $(TARGET)
	./$(TARGET)

debug: $(TARGET)
	valgrind ./$(TARGET)
