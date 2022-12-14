PREFIX = /usr/local/bin/
INC=-Isrc -Isrc/commands -Isrc/index -Isrc/shared -Isrc/transform -Isrc/search
VPATH=src:src/commands:src/index:src/shared:src/transform:src/search

SRCS =  \
	alignment.cpp \
	assemble.cpp \
	consensible.cpp \
	consensus.cpp \
	filenames.cpp \
	index.cpp \
	index_reader.cpp \
	index_structs.cpp \
	index_writer.cpp \
	local_alignment.cpp \
	match.cpp \
	match_query.cpp \
	parameters.cpp \
	query_binary.cpp \
	query_structs.cpp \
	read.cpp \
	result.cpp \
	shared_functions.cpp \
	shared_structs.cpp \
	target.cpp \
	timer.cpp \
	transform.cpp \
	transform_bwt.cpp \
	transform_structs.cpp \
	transform_binary.cpp


# C++ compiler
CXX = g++
# C++ flags; passed to compiler
CXXFLAGS = -std=c++11
# Linker flags; passed to compiler
LDFLAGS = -std=c++11
# Dependency flags; passed to compiler
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
# Objects directory
OBJDIR = .o
$(shell mkdir -p $(OBJDIR) >/dev/null)
# Dependencies directory
DEPDIR = .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
# Derive objects from sources
OBJS = $(patsubst %,$(OBJDIR)/%.o,$(basename $(SRCS)))
# Derive dependencies from sources
DEPS = $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))

# Generic link executable
LINK.o = $(CXX) $(LDFLAGS) $(DBG) -o $@
# Generic compile object
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(DBG) $(INC) -c -o $@
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

all: consensible
	@echo
	@echo 'Compile successful. Type "sudo make install" to complete intall.'

.PHONY: install
install:
	@mv -f consensible $(PREFIX)
	@make clean
	@echo
	@echo 'Install successful. Type "consensible -h" to see usage.'

.PHONY: clean
clean:
	@$(RM) -r $(OBJDIR) $(DEPDIR)

consensible: $(OBJS)
	$(LINK.o) $^

$(OBJDIR)/%.o : %.cpp
$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d
	$(COMPILE.cc) $<
	$(POSTCOMPILE)

# Create dependency rule
.PRECIOUS = $(DEPDIR)/%.d
$(DEPDIR)/%.d: ;

# Include dependencies; The '-' ensures no errors
-include $(DEPS)
