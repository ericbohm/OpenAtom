include ../config.mk

C_FILES=$(wildcard *.C)
CI_FILES=$(wildcard *.ci)
MODULES=$(CI_FILES:.ci=.decl.h)
OBJS=$(addprefix $(BUILDDIR)/, $(C_FILES:.C=.o))

default: all
all: $(BUILDDIR) $(OBJS)

$(BUILDDIR):
	mkdir $(BUILDDIR)

make_ci: $(MODULES)

%.decl.h: %.ci
	$(CHARMC) $<

$(BUILDDIR)/%.o: %.C %.h $(MODULES)
	$(CHARMC) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.C $(MODULES)
	$(CHARMC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf *.decl.h *.def.h
